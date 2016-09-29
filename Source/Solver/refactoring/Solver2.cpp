//
// Created by Dmitri Bagaev on 22.09.16.
//

#include <Source/Solver/refactoring/solver_petsc/SolverPETSc.h>
#include <inmost_sparse.h>
#include "Solver2.h"

namespace INMOST {

    int *Solver2::argc = NULL;
    char ***Solver2::argv = NULL;
    const char *Solver2::database = NULL;
    bool Solver2::is_initialized = false;
    bool Solver2::is_finalized = false;

    Solver2::Solver2(std::string solverName, std::string prefix, INMOST_MPI_Comm _comm) {
        this->solver = SolverFactory::getSolver(solverName);
        this->prefix = prefix;

        solver->SetCommunicator(_comm);
        std::string solverDatabasePath = Solver2::parseDatabase(solverName);
        solver->Initialize(argc, argv, solverDatabasePath.c_str(), prefix);
    }

    Solver2::Solver2(const Solver2 &other) {
        this->solver = SolverFactory::copySolver(other.solver);
        this->prefix = other.prefix;

        solver->SetCommunicator(other.solver->GetCommunicator());
        std::string solverDatabasePath = Solver2::parseDatabase(solver->SolverName());
        solver->Initialize(argc, argv, solverDatabasePath.c_str(), prefix);
    }

    Solver2& Solver2::operator=(const Solver2& other) {
        if( this != &other ) {
            this->solver->SetCommunicator(other.solver->GetCommunicator());
            this->prefix = other.prefix;
            this->solver->Assign(other.solver);
        }
        return *this;
    }

    void Solver2::Initialize(int *argc, char ***argv, const char *database) {
        Solver2::argc = argc;
        Solver2::argv = argv;
        Solver2::database = database;
        Solver2::is_initialized = true;
        Solver2::is_finalized = false;
        Sparse::CreateRowEntryType();
        //Register all available solvers
#if defined(USE_SOLVER_PETSC)
        SolverFactory::registerSolver<SolverPETSc>("petsc");
#endif
    }

    void Solver2::SetMatrix(Sparse::Matrix & A, bool ModifiedPattern, bool OldPreconditioner) {
        solver->SetMatrix(A, ModifiedPattern, OldPreconditioner);
    }

    bool Solver2::Solve(INMOST::Sparse::Vector & RHS, INMOST::Sparse::Vector & SOL) {
        if( !solver->isMatrixSet()) throw MatrixNotSetInSolver;
        if( RHS.GetCommunicator() != solver->GetCommunicator() || SOL.GetCommunicator() != solver->GetCommunicator()) throw DifferentCommunicatorInSolver;
        INMOST_DATA_ENUM_TYPE vbeg,vend;
        RHS.GetInterval(vbeg,vend);
        if( RHS.Size() != SOL.Size() )
        {
            if( SOL.Size() == 0 )
            {
                SOL.SetInterval(vbeg,vend);
                for(Sparse::Vector::iterator ri = SOL.Begin(); ri != SOL.End(); ++ri) *ri = 0.0;
            }
            else throw InconsistentSizesInSolver;
        }
        return solver->Solve(RHS, SOL);
    }

    void Solver2::Finalize() {
        Sparse::DestroyRowEntryType();
        Solver2::is_finalized = true;
        Solver2::is_initialized = false;
    }

    bool Solver2::isInitialized() {
        return is_initialized;
    }

    bool Solver2::isFinalized() {
        return is_finalized;
    }

    const INMOST_DATA_ENUM_TYPE Solver2::Iterations() const {
        return solver->Iterations();
    }

    const INMOST_DATA_REAL_TYPE Solver2::Residual() const {
        return solver->Residual();
    }

    const std::string Solver2::ReturnReason() const {
        return solver->ReturnReason();
    }

    std::string Solver2::SolverName() const {
        return solver->SolverName();
    }

    std::string Solver2::SolverPrefix() const {
        return prefix;
    }

    Solver2::~Solver2() {
        solver->Finalize();
        delete solver;
    }

    std::string Solver2::parseDatabase(std::string solverName) {
        const char *name = solverName.c_str();
        if( database != NULL ) {
            FILE * f = fopen(database, "r");
            if (f != NULL) {
                char str[4096];
                while( !feof(f) && fgets(str,4096,f)) {
                    int k = 0, l;
                    for(k = 0; k < (int)strlen(str); ++k) {
                        if( str[k] == ':' ) break;
                    }
                    if( k == strlen(str) ) continue; //invalid line
                    for(l = 0; l < k; ++l) str[l] = tolower(str[l]);
                    l = (int)strlen(str)-1; // Right-trim string
                    while(l > 0 && isspace(str[l]) ) --l;
                    str[l+1] = 0;
                    l = k+1;
                    while(l < (int)strlen(str) && isspace(str[l]) ) ++l;
                    if( l == strlen(str) ) continue; //skip empty entry
                    if( !strncmp(str, name, k) ) {
                        return std::string(str+l);
                    }
                }
                fclose(f);
            }
        }
        return std::string("");
    }

}
