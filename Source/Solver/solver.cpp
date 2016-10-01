#include <inmost.h>

#include "solver_inner/solver_ilu2/SolverILU2.h"
#include "solver_inner/solver_ddpqiluc2/SolverDDPQILUC2.h"
#include "solver_inner/solver_mptiluc/SolverMPTILUC.h"
#include "solver_inner/solver_mptilu2/SolverMPTILU2.h"

#if defined(USE_SOLVER_PETSC)
#include "solver_petsc/SolverPETSc.h"
#endif

#if defined(USE_SOLVER_TRILINOS)
#include "solver_trilinos/SolverTrilinos.h"
#include "solver_trilinos/solver_aztec/SolverTrilinosAztec.h"
#include "solver_trilinos/solver_belos/SolverTrilinosBelos.h"
#include "solver_trilinos/solver_ml/SolverTrilinosML.h"
#include "solver_trilinos/solver_ifpack/SolverTrilinosIfpack.h"
#endif

#if defined(USE_SOLVER_ANI)
#include "solver_ani/SolverANI.h"
#endif

namespace INMOST {

    int *Solver::argc = NULL;
    char ***Solver::argv = NULL;
    const char *Solver::database = NULL;
    bool Solver::is_initialized = false;
    bool Solver::is_finalized = false;

    Solver::Solver(std::string solverName, std::string prefix, INMOST_MPI_Comm _comm) {
        this->solver = SolverMaster::getSolver(solverName);
        this->prefix = prefix;

        solver->SetCommunicator(_comm);
        std::string solverDatabasePath = Solver::parseDatabase(solverName);
        solver->Initialize(argc, argv, solverDatabasePath.c_str(), prefix);
    }

    Solver::Solver(const Solver &other) {
        this->solver = SolverMaster::copySolver(other.solver);
        this->prefix = other.prefix;

        solver->SetCommunicator(other.solver->GetCommunicator());
        std::string solverDatabasePath = Solver::parseDatabase(solver->SolverName());
        solver->Initialize(argc, argv, solverDatabasePath.c_str(), prefix);
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
        Solver::database = database;
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
#if defined(USE_SOLVER_TRILINOS)
        SolverMaster::registerSolver<SolverTrilinosAztec>("trilinos_aztec");
        SolverMaster::registerSolver<SolverTrilinosBelos>("trilinos_belos");
        SolverMaster::registerSolver<SolverTrilinosML>("trilinos_ml");
        SolverMaster::registerSolver<SolverTrilinosIfpack>("trilinos_ifpack");
#endif
#if defined(USE_SOLVER_ANI)
        SolverMaster::registerSolver<SolverANI>("ani");
#endif
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

    INMOST_DATA_REAL_TYPE Solver::GetParameterReal(std::string property) const {
        return solver->GetParameterReal(property);
    }

    INMOST_DATA_ENUM_TYPE Solver::GetParameterEnum(std::string property) const {
        return solver->GetParameterEnum(property);
    }

    void Solver::SetParameterReal(std::string property, INMOST_DATA_REAL_TYPE value) {
        solver->SetParameterReal(property, value);
    }

    void Solver::SetParameterEnum(std::string property, INMOST_DATA_ENUM_TYPE value) {
        solver->SetParameterEnum(property, value);
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

    Solver::~Solver() {
        solver->Finalize();
        delete solver;
    }

    std::string Solver::parseDatabase(std::string solverName) {
        const char *name = solverName.c_str();
        if (database != NULL) {
            FILE *f = fopen(database, "r");
            if (f != NULL) {
                char str[4096];
                while (!feof(f) && fgets(str, 4096, f)) {
                    int k = 0, l;
                    for (k = 0; k < (int) strlen(str); ++k) {
                        if (str[k] == ':') break;
                    }
                    if (k == strlen(str)) continue; //invalid line
                    for (l = 0; l < k; ++l) str[l] = tolower(str[l]);
                    l = (int) strlen(str) - 1; // Right-trim string
                    while (l > 0 && isspace(str[l])) --l;
                    str[l + 1] = 0;
                    l = k + 1;
                    while (l < (int) strlen(str) && isspace(str[l])) ++l;
                    if (l == strlen(str)) continue; //skip empty entry
                    if (!strncmp(str, name, k)) {
                        return std::string(str + l);
                    }
                }
                fclose(f);
            }
        }
        return std::string("");
    }

}