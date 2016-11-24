#include <inmost.h>

#include "SolverMaster.h"
#include "solver_inner/solver_ilu2/SolverILU2.h"
#include "solver_inner/solver_ddpqiluc2/SolverDDPQILUC2.h"
#include "solver_inner/solver_mptiluc/SolverMPTILUC.h"
#include "solver_inner/solver_mptilu2/SolverMPTILU2.h"

#if defined(USE_SOLVER_PETSC)

#include "solver_petsc/SolverPETSc.h"

#endif

#if defined(USE_SOLVER_TRILINOS) && defined(USE_MPI)

#include "solver_trilinos/SolverTrilinos.h"
#include "solver_trilinos/solver_aztec/SolverTrilinosAztec.h"
#include "solver_trilinos/solver_belos/SolverTrilinosBelos.h"
#include "solver_trilinos/solver_ml/SolverTrilinosML.h"
#include "solver_trilinos/solver_ifpack/SolverTrilinosIfpack.h"

#endif

#if defined(USE_SOLVER_ANI)

#include "solver_ani/SolverANI.h"

#endif

#if defined(USE_SOLVER_SUPERLU)

#include "solver_superlu/SolverSUPERLU.h"

#endif

#if defined(HAVE_SOLVER_K3BIILU2)

#include "solver_k3biilu2/SolverK3BIILU2.h"

#endif

#if defined(HAVE_SOLVER_FCBIILU2)

#include "solver_fcbiilu2/SolverFCBIILU2.h"

#endif

namespace INMOST {

    int *Solver::argc = NULL;
    char ***Solver::argv = NULL;
    bool Solver::is_initialized = false;
    bool Solver::is_finalized = false;
    std::vector<SolverParameters> Solver::parameters = std::vector<SolverParameters>();

    Solver::Solver(std::string solverName, std::string prefix, INMOST_MPI_Comm _comm) {
        std::string lowerName = string_to_lower(solverName);
        this->solver = SolverMaster::getSolver(lowerName);
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

    Solver::Solver(const Solver &other) {
        this->solver = SolverMaster::copySolver(other.solver);
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
            this->prefix = other.prefix;
            this->solver->Assign(other.solver);
        }
        return *this;
    }

    std::string Solver::SolverName() const {
        return solver->SolverName();
    }

    std::string Solver::SolverPrefix() const {
        return prefix;
    }

    void Solver::Initialize(int *argc, char ***argv, const char *database) {
        Solver::argc = argc;
        Solver::argv = argv;
        Solver::is_initialized = true;
        Solver::is_finalized = false;
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
        //Register all available solvers
        SolverMaster::registerSolver<SolverILU2>("inner_ilu2");
        SolverMaster::registerSolver<SolverDDPQILUC2>("inner_ddpqiluc2");
        SolverMaster::registerSolver<SolverMPTILUC>("inner_mptiluc");
        SolverMaster::registerSolver<SolverMPTILU2>("inner_mptilu2");
#if defined(USE_SOLVER_PETSC)
        SolverMaster::registerSolver<SolverPETSc>("petsc");
#endif
#if defined(USE_SOLVER_TRILINOS) && defined(USE_MPI)
        SolverMaster::registerSolver<SolverTrilinosAztec>("trilinos_aztec");
        SolverMaster::registerSolver<SolverTrilinosBelos>("trilinos_belos");
        SolverMaster::registerSolver<SolverTrilinosML>("trilinos_ml");
        SolverMaster::registerSolver<SolverTrilinosIfpack>("trilinos_ifpack");
#endif
#if defined(USE_SOLVER_ANI)
        SolverMaster::registerSolver<SolverANI>("ani");
#endif
#if defined(USE_SOLVER_SUPERLU)
        SolverMaster::registerSolver<SolverSUPERLU>("superlu");
#endif
#if defined(HAVE_SOLVER_K3BIILU2)
        SolverMaster::registerSolver<SolverK3BIILU2>("k3biilu2");
#endif
#if defined(HAVE_SOLVER_FCBIILU2)
        SolverMaster::registerSolver<SolverFCBIILU2>("fcbiilu2");
#endif
        Solver::parseXMLDatabase(database);

        //Debug
//        for (auto p = parameters.begin(); p < parameters.end(); p++) {
//            std::cout << "============================================================================" << std::endl;
//            std::cout << (*p).solverName << ":" << (*p).solverPrefix << ":" << (*p).internalFile << std::endl;
//            for (auto pp = (*p).parameters.begin(); pp < (*p).parameters.end(); pp++) {
//                std::cout << (*pp).first << " = " << (*pp).second << std::endl;
//            }
//            std::cout << "============================================================================" << std::endl;
//        }

        Sparse::CreateRowEntryType();
    }

    void Solver::Finalize() {
        Sparse::DestroyRowEntryType();
#if defined(USE_MPI)
        {
            int flag = 0;
            MPI_Finalized(&flag);
            if (!flag) {
                MPI_Finalize();
            }
        }
#endif
        Solver::is_finalized = true;
        Solver::is_initialized = false;
    }

    bool Solver::isInitialized() {
        return is_initialized;
    }

    bool Solver::isFinalized() {
        return is_finalized;
    }

    void Solver::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) {
        solver->SetMatrix(A, ModifiedPattern, OldPreconditioner);
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
        return solver->Solve(RHS, SOL);
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

    const INMOST_DATA_ENUM_TYPE Solver::Iterations() const {
        return solver->Iterations();
    }

    const INMOST_DATA_REAL_TYPE Solver::Residual() const {
        return solver->Residual();
    }

    const std::string Solver::ReturnReason() const {
        return solver->ReturnReason();
    }

    INMOST_DATA_REAL_TYPE Solver::Condest(INMOST_DATA_REAL_TYPE tol, unsigned int maxits) {
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
        if (xml_database == NULL) return;

        std::ifstream input;
        input.open(xml_database);

        if (input.fail()) {
            std::cout << __FILE__ << ": XML database file not found " << std::endl;
            return;
        }

        XMLReader reader(std::string(xml_database), input);

        try {
            XMLReader::XMLTree root = reader.ReadXML();

            if (root.tag.name != "SolverParameters") {
                std::cout << __FILE__ << ": Bad XML database file" << std::endl;
                return;
            }

            for (xml_reader_tree_iterator_t solver = root.children.begin(); solver < root.children.end(); solver++) {
                std::string solverName = string_to_lower((*solver).tag.name);
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

}