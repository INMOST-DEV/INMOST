#ifndef INMOST_SOLVER_MASTER
#define INMOST_SOLVER_MASTER

#include <inmost_solver.h>
#include <inmost_solver_interface.h>

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

    class SolverMaster {
    private:
        static SolverInterface *getSolver(std::string name);

        static std::vector<std::string> getAvailableSolvers();

        static bool isSolverAvailable(std::string name);

        friend Solver::Solver(std::string solverName, std::string prefix, INMOST_MPI_Comm _comm);

        friend Solver::Solver(const Solver &other);

        friend void Solver::Initialize(int *argc, char ***argv, const char *database);

        friend bool Solver::isSolverAvailable(std::string name);

        friend std::vector<std::string> Solver::getAvailableSolvers();
    };
}

#endif //INMOST_SOLVER_MASTER
