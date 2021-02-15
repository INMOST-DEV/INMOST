#include <inmost.h>
#include "SolverMaster.h"
#include "../Misc/utils.h"

#if defined(USE_SOLVER)

namespace INMOST {

    int  *Solver::argc          = NULL;
    char ***Solver::argv        = NULL;
    bool Solver::is_initialized = false;
    bool Solver::is_finalized   = false;
    std::vector<SolverParameters> Solver::parameters = std::vector<SolverParameters>();

    void Solver::SetParameterEnum(std::string name, INMOST_DATA_ENUM_TYPE value) { SetParameter(name, to_string(value)); }

    void Solver::SetParameterReal(std::string name, INMOST_DATA_REAL_TYPE value) { SetParameter(name, to_string(value)); }

    Solver::Solver(std::string solverName, std::string prefix, INMOST_MPI_Comm _comm) :
            versbosity(SolverVerbosityLevel0), preconditioner_time(0), iterations_time(0), is_solved(false) {
        std::string lowerName = string_to_lower(solverName);
        this->solver = SolverMaster::getSolver(solverName);
        this->prefix = string_to_lower(prefix);
        solver->SetCommunicator(_comm);
        //TODO find easiest way
        bool parametersFound = false;
        if (Solver::parameters.size() > 0) {
            for (solver_parameters_iterator_t parameters = Solver::parameters.end() - 1; parameters >= Solver::parameters.begin(); parameters--) {
                if ((*parameters).solverName == lowerName && (*parameters).solverPrefix == (this->prefix)) {
                    solver->Setup(argc, argv, *parameters);
                    parametersFound = true;
                    break;
                }
            }
        }
        if (!parametersFound) {
            SolverParameters emptyParameters(lowerName, this->prefix, "");
            solver->Setup(argc, argv, emptyParameters);
        }

    }

    Solver::Solver(const Solver &other) :
            versbosity(other.versbosity), preconditioner_time(other.preconditioner_time), iterations_time(other.iterations_time), is_solved(other.is_solved) {
        this->solver = SolverMaster::getSolver(other.solver->SolverName());
        this->solver->Copy(other.solver);
        this->prefix = other.prefix;
        solver->SetCommunicator(other.solver->GetCommunicator());
        //TODO find easiest way
        bool parametersFound = false;
        if (Solver::parameters.size() > 0) {
            for (solver_parameters_iterator_t parameters = Solver::parameters.end() - 1; parameters >= Solver::parameters.begin(); parameters--) {
                if ((*parameters).solverName == other.solver->SolverName() && (*parameters).solverPrefix == (this->prefix)) {
                    solver->Setup(argc, argv, *parameters);
                    parametersFound = true;
                    break;
                }
            }
        }
        if (!parametersFound) {
            SolverParameters emptyParameters(other.solver->SolverName(), this->prefix, "");
            solver->Setup(argc, argv, emptyParameters);
        }
    }

    Solver &Solver::operator=(const Solver &other) {
        if (this != &other) {
            this->solver->SetCommunicator(other.solver->GetCommunicator());
            this->solver->Assign(other.solver);

            this->prefix              = other.prefix;
            this->versbosity          = other.versbosity;
            this->preconditioner_time = other.preconditioner_time;
            this->iterations_time     = other.iterations_time;
            this->is_solved           = other.is_solved;
        }
        return *this;
    }

    std::string Solver::SolverName() const {
        return solver->SolverName();
    }

    std::string Solver::SolverPrefix() const {
        return prefix;
    }

    SolverVerbosityLevel Solver::VerbosityLevel() const {
        return versbosity;
    }

    void Solver::SetVerbosityLevel(INMOST::SolverVerbosityLevel level) {
        versbosity = level;
    }

    void Solver::Initialize(int *argc, char ***argv, const char *database) {
        Solver::argc           = argc;
        Solver::argv           = argv;
        Solver::is_initialized = true;
        Solver::is_finalized   = false;
#if defined(USE_MPI)
        {
            int flag = 0;
            int ierr = 0;
            MPI_Initialized(&flag);
            if (!flag) {
                ierr = MPI_Init(argc, argv);
                if (ierr != MPI_SUCCESS) {
                    std::cout << __FILE__ << ":" << __LINE__ << "problem in MPI_Init" << std::endl;
                }
            }
        }
#endif
        Solver::parseXMLDatabase(database);
        Sparse::CreateRowEntryType();
    }

    void Solver::Finalize() {
        Sparse::DestroyRowEntryType();
#if defined(USE_MPI)
        {
            int flag = 0;
            MPI_Finalized(&flag);
            if (!flag)
                MPI_Finalize();
        }
#endif
        Solver::is_finalized   = true;
        Solver::is_initialized = false;
    }

    bool Solver::isInitialized() {
        return is_initialized;
    }

    bool Solver::isFinalized() {
        return is_finalized;
    }

    void Solver::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) {
        double preconditioner_timer_start = Timer();
        solver->SetMatrix(A, ModifiedPattern, OldPreconditioner);
        //INMOST::MPIBarrier();
        preconditioner_time = Timer() - preconditioner_timer_start;
    }

    bool Solver::Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL) {
        if (!solver->isMatrixSet()) throw MatrixNotSetInSolver;
        if (RHS.GetCommunicator() != solver->GetCommunicator() ||
            SOL.GetCommunicator() != solver->GetCommunicator())
            throw DifferentCommunicatorInSolver;
        INMOST_DATA_ENUM_TYPE vbeg, vend;
        RHS.GetInterval(vbeg, vend);
        if (RHS.Size() != SOL.Size()) {
            if (SOL.Size() == 0) {
                SOL.SetInterval(vbeg, vend);
                for (Sparse::Vector::iterator ri = SOL.Begin(); ri != SOL.End(); ++ri) *ri = 0.0;
            } else throw InconsistentSizesInSolver;
        }

        double iterations_timer_start = Timer();

        is_solved = solver->Solve(RHS, SOL);
        //INMOST::MPIBarrier();
        iterations_time = Timer() - iterations_timer_start;

        if (versbosity > SolverVerbosityLevel1) {
#if defined(USE_MPI)
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if (rank == 0) {
                std::cout << SolutionMetadataLine("\t") << std::endl;
            }
#endif
        }

        return is_solved;
    }

    bool Solver::Clear() {
        return solver->Clear();
    }

    std::string Solver::GetParameter(std::string name) const {
        return solver->GetParameter(name);
    }

    void Solver::SetParameter(std::string name, std::string value) {
        return solver->SetParameter(name, value);
    }

    INMOST_DATA_ENUM_TYPE Solver::Iterations() const {
        return solver->Iterations();
    }

    INMOST_DATA_REAL_TYPE Solver::Residual() const {
        return solver->Residual();
    }

    const std::string Solver::ReturnReason() const {
        return solver->ReturnReason();
    }

    double Solver::PreconditionerTime() const {
        return preconditioner_time;
    }

    double Solver::IterationsTime() const {
        return iterations_time;
    }

    bool Solver::IsSolved() const {
        return is_solved;
    }

    const std::string Solver::SolutionMetadataLine(const std::string &delimiter) const {
        return this->SolverName() + delimiter +
               this->SolverPrefix() + delimiter +
               INMOST::to_string(this->IsSolved()) + delimiter +
               INMOST::to_string(this->Iterations()) + delimiter +
               INMOST::to_string(this->Residual()) + delimiter +
               INMOST::to_string(this->PreconditionerTime()) + delimiter +
               INMOST::to_string(this->IterationsTime()) + delimiter +
               INMOST::to_string(this->PreconditionerTime() + this->IterationsTime());
    }

    INMOST_DATA_REAL_TYPE Solver::Condest(INMOST_DATA_REAL_TYPE tol, INMOST_DATA_ENUM_TYPE maxits) {
        if (!solver->isMatrixSet()) throw MatrixNotSetInSolver;
        return solver->Condest(tol, maxits);
    }

    bool Solver::isSolverAvailable(std::string name) {
        return SolverMaster::isSolverAvailable(name);
    }

    std::vector<std::string> Solver::getAvailableSolvers() {
        return SolverMaster::getAvailableSolvers();
    }

    Solver::~Solver() {
        solver->Finalize();
        delete solver;
    }

    void Solver::parseXMLDatabase(const char *xml_database) {
        if (xml_database == NULL || xml_database[0] == '\0') return;

        std::ifstream input;
        input.open(xml_database);

        if (input.fail()) {
            std::cout << __FILE__ << ":" << __LINE__ << ": XML database file " << xml_database << " not found." << std::endl;
            return;
        }

        XMLReader reader(std::string(xml_database), input);

        try {
            XMLReader::XMLTree root = reader.ReadXML();

            if (root.tag.name != "SolverParameters") {
                std::cout << __FILE__ << ":" << __LINE__ << ": Bad XML database file " << xml_database << "!" << std::endl;
                return;
            }

            for (xml_reader_tree_iterator_t solver = root.children.begin(); solver < root.children.end(); solver++) {
                std::string solverName   = string_to_lower((*solver).tag.name);
                std::string internalFile = "";

                if ((*solver).tag.attributes.size() != 0) {
                    for (xml_reader_attrib_iterator_t attr = (*solver).tag.attributes.begin(); attr < (*solver).tag.attributes.end(); attr++) {
                        if ((*attr).name == "File" || (*attr).name == "file" || (*attr).name == "FILE") {
                            internalFile = (*attr).value;
                        }
                    }
                }

                if ((*solver).children.size() == 0) {
                    //Internal solver
                    parameters.push_back(SolverParameters(solverName, "", internalFile));
                } else {
                    //Inner solver
                    for (xml_reader_tree_iterator_t prefix = (*solver).children.begin(); prefix < (*solver).children.end(); prefix++) {
                        internalFile = "";
                        std::string solverPrefix = string_to_lower((*prefix).tag.name);

                        if ((*prefix).tag.attributes.size() != 0) {
                            for (xml_reader_attrib_iterator_t attr = (*prefix).tag.attributes.begin(); attr < (*prefix).tag.attributes.end(); attr++) {
                                if ((*attr).name == "File" || (*attr).name == "file" || (*attr).name == "FILE") {
                                    internalFile = (*attr).value;
                                }
                            }
                        }

                        SolverParameters prefix_p = SolverParameters(solverName, solverPrefix, internalFile);

                        for (xml_reader_tree_iterator_t p = (*prefix).children.begin(); p < (*prefix).children.end(); p++) {

                            if ((*p).tag.attributes.size() == 1) {
                                if ((*p).tag.attributes[0].name == "value" || (*p).tag.attributes[0].name == "Value") {
                                    prefix_p.parameters.push_back(std::make_pair((*p).tag.name, (*p).tag.attributes[0].value));
                                }
                            }
                        }

                        parameters.push_back(prefix_p);

                    }

                }

            }

        } catch (...) {
            std::cout << __FILE__ << ": Error while parsing xml database file" << std::endl;
        }

    }

	INMOST_DATA_REAL_TYPE Sparse::Matrix::Residual(Sparse::Vector &RHS, Sparse::Vector &SOL)
	{
		Sparse::Vector R(RHS);
		INMOST_DATA_ENUM_TYPE mbeg, mend, k;
		if( isParallel() )
		{
			std::cout << __FILE__ << ":" << __LINE__ << " matrix already in parallel state " << std::endl;
			MatVec(1.0, SOL, -1.0, R);
			mbeg = GetFirstIndex();
			mend = GetLastIndex();
#if defined(USE_MPI)
			int rank, size;
			MPI_Comm_rank(GetCommunicator(),&rank);
			MPI_Comm_size(GetCommunicator(),&size);
			if( rank == size-1 )
				MPI_Send(&mend,1,INMOST_MPI_DATA_ENUM_TYPE,rank-1,rank-1,
						 GetCommunicator());//,MPI_STATUS_IGNORE);
			else if( rank > 0 )
				MPI_Sendrecv(&mbeg,1,INMOST_MPI_DATA_ENUM_TYPE,rank-1,rank-1,
							 &mend,1,INMOST_MPI_DATA_ENUM_TYPE,rank+1,rank,
							 GetCommunicator(),MPI_STATUS_IGNORE);
			else
				MPI_Recv(&mend,1,INMOST_MPI_DATA_ENUM_TYPE,rank+1,rank,
						 GetCommunicator(),MPI_STATUS_IGNORE);
			std::cout << "on " << rank << " " << mbeg << ":" << mend << std::endl;
#endif
		}
		else
		{
			Solver::OrderInfo info;
			info.PrepareMatrix(*this, 0);
			info.PrepareVector(SOL);
			info.GetLocalRegion(info.GetRank(), mbeg, mend);
			info.Update(SOL);
			MatVec(1.0, SOL, -1.0, R);
			info.RestoreVector(SOL);
			info.RestoreMatrix(*this);

		}
		INMOST_DATA_REAL_TYPE aresid = 0, bresid = 0;
		for (k = mbeg; k < mend; ++k)
		{
			aresid += R[k] * R[k];
			bresid += RHS[k] * RHS[k];
		}
		double temp[2] = {aresid, bresid}, recv[2] = {aresid, bresid};
#if defined(USE_MPI)
		MPI_Allreduce(temp, recv, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
		return sqrt(recv[0]/(recv[1] + 1.0e-100));
	}
}

#endif //USE_SOLVER
