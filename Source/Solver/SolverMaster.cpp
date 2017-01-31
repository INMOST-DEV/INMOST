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
	
	//initialization of backward-compatibility constants
	const Solver::Type Solver::INNER_ILU2 = "inner_ilu2";
	const Solver::Type Solver::INNER_DDPQILUC = "inner_ddpqiluc2";
	const Solver::Type Solver::INNER_MPTILUC = "inner_mptiluc";
	const Solver::Type Solver::INNER_MPTILU2 = "inner_mptilu2";
	const Solver::Type Solver::Trilinos_Aztec = "trilinos_aztec";
	const Solver::Type Solver::Trilinos_Belos = "trilinos_belos";
	const Solver::Type Solver::Trilinos_ML = "trilinos_ml";
	const Solver::Type Solver::Trilinos_Ifpack = "trilinos_ifpack";
	const Solver::Type Solver::PETSc = "petsc";
	const Solver::Type Solver::ANI = "inner_ilu2";
	const Solver::Type Solver::FCBIILU2 = "fcbiilu2";
	const Solver::Type Solver::K3BIILU2 = "k3biilu2";
	const Solver::Type Solver::SUPERLU = "superlu";

    SolverInterface *SolverMaster::getSolver(std::string name) {
        if (name == "inner_ilu2") return new SolverILU2();
        if (name == "inner_ddpqiluc2") return new SolverDDPQILUC2();
        if (name == "inner_mptiluc") return new SolverMPTILUC();
        if (name == "inner_mptilu2") return new SolverMPTILU2();
#if defined(USE_SOLVER_PETSC)
        if (name == "petsc") return new SolverPETSc();
#endif
#if defined(USE_SOLVER_TRILINOS) && defined(USE_MPI)
        if (name == "trilinos_aztec") return new SolverTrilinosAztec();
        if (name == "trilinos_belos") return new SolverTrilinosBelos();
        if (name == "trilinos_ml") return new SolverTrilinosML();
        if (name == "trilinos_ifpack") return new SolverTrilinosIfpack();
#endif
#if defined(USE_SOLVER_ANI)
        if (name == "ani") return new SolverANI();
#endif
#if defined(USE_SOLVER_SUPERLU)
        if (name == "superlu") return new SolverSUPERLU();
#endif
#if defined(HAVE_SOLVER_K3BIILU2)
        if (name == "k3biilu2") return new SolverK3BIILU2();
#endif
#if defined(HAVE_SOLVER_FCBIILU2)
        if (name == "fcbiilu2") return new SolverFCBIILU2();
#endif
        throw INMOST::SolverNotFound;
    }

    std::vector<std::string> SolverMaster::getAvailableSolvers() {
        std::vector<std::string> s;
        s.push_back("inner_ilu2");
        s.push_back("inner_ddpqiluc2");
        s.push_back("inner_mptiluc");
        s.push_back("inner_mptilu2");
#if defined(USE_SOLVER_PETSC)
        s.push_back("petsc");
#endif
#if defined(USE_SOLVER_TRILINOS) && defined(USE_MPI)
        s.push_back("trilinos_aztec");
        s.push_back("trilinos_belos");
        s.push_back("trilinos_ml");
        s.push_back("trilinos_ifpack");
#endif
#if defined(USE_SOLVER_ANI)
        s.push_back("ani");
#endif
#if defined(USE_SOLVER_SUPERLU)
        s.push_back("superlu");
#endif
#if defined(HAVE_SOLVER_K3BIILU2)
        s.push_back("k3biilu2");
#endif
#if defined(HAVE_SOLVER_FCBIILU2)
        s.push_back("fcbiilu2");
#endif
        return s;
    }

    bool SolverMaster::isSolverAvailable(std::string name) {
        try {
            SolverInterface *s = SolverMaster::getSolver(name);
            delete s;
            return true;
        } catch (...) {
            return false;
        }
    }

}

