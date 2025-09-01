#include "SolverMaster.h"
#if defined(USE_SOLVER)
#include "solver_inner/solver_ilu2/SolverILU2.h"
#include "solver_inner/solver_ddpqiluc2/SolverDDPQILUC2.h"
#include "solver_inner/solver_mptiluc/SolverMPTILUC.h"
#include "solver_inner/solver_mlmptiluc/SolverMLMPTILUC.h"
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

#if defined(USE_SOLVER_HYPRE)

#include "solver_hypre/SolverHYPRE.h"

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
	const Solver::Type Solver::INNER_MLMPTILUC = "inner_mlmptiluc";
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
	const Solver::Type Solver::HYPRE_BOOMERAMG = "hypre_boomeramg";
    const Solver::Type Solver::HYPRE_PARASAILS = "hypre_parasails";
    const Solver::Type Solver::HYPRE_PILUT = "hypre_pilut";
    const Solver::Type Solver::HYPRE_AMS = "hypre_ams";
    const Solver::Type Solver::HYPRE_EUCLID = "hypre_euclid";

    SolverInterface *SolverMaster::getSolver(std::string name)
	{
        if (name == Solver::INNER_ILU2) return new SolverILU2();
        if (name == Solver::INNER_DDPQILUC) return new SolverDDPQILUC2();
        if (name == Solver::INNER_MPTILUC) return new SolverMPTILUC();
		if (name == Solver::INNER_MLMPTILUC) return new SolverMLMPTILUC();
        if (name == Solver::INNER_MPTILU2) return new SolverMPTILU2();
        
#if defined(USE_SOLVER_PETSC)
        if (name == Solver::PETSc) return new SolverPETSc();
#endif
        
#if defined(USE_SOLVER_TRILINOS) && defined(USE_MPI)
        if (name == Solver::Trilinos_Aztec) return new SolverTrilinosAztec();
        if (name == Solver::Trilinos_Belos) return new SolverTrilinosBelos();
        if (name == Solver::Trilinos_ML) return new SolverTrilinosML();
        if (name == Solver::Trilinos_Ifpack) return new SolverTrilinosIfpack();
#endif
        
#if defined(USE_SOLVER_ANI)
        if (name == Solver::ANI) return new SolverANI();
#endif
        
#if defined(USE_SOLVER_SUPERLU)
        if (name == Solver::SUPERLU) return new SolverSUPERLU();
#endif

#if defined(USE_SOLVER_HYPRE)
        if (name == Solver::HYPRE_BOOMERAMG) return new SolverHYPRE(BoomerAMG);
        if (name == Solver::HYPRE_PARASAILS) return new SolverHYPRE(ParaSails);
        if (name == Solver::HYPRE_PILUT) return new SolverHYPRE(PILUT);
        if (name == Solver::HYPRE_AMS) return new SolverHYPRE(AMS);
        if (name == Solver::HYPRE_EUCLID) return new SolverHYPRE(Euclid);
#endif
        
#if defined(HAVE_SOLVER_K3BIILU2)
        if (name == Solver::K3BIILU2) return new SolverK3BIILU2();
#endif
#if defined(HAVE_SOLVER_FCBIILU2)
        if (name == Solver::FCBIILU2) return new SolverFCBIILU2();
#endif
        throw INMOST::SolverNotFound;
    }

    std::vector<std::string> SolverMaster::getAvailableSolvers()
	{
        std::vector<std::string> s;
        s.push_back(Solver::INNER_ILU2);
        s.push_back(Solver::INNER_DDPQILUC);
        s.push_back(Solver::INNER_MLMPTILUC);
		s.push_back(Solver::INNER_MLMPTILUC);
        s.push_back(Solver::INNER_MPTILU2);
#if defined(USE_SOLVER_PETSC)
        s.push_back(Solver::PETSc);
#endif
        
#if defined(USE_SOLVER_TRILINOS) && defined(USE_MPI)
        s.push_back(Solver::Trilinos_Aztec);
        s.push_back(Solver::Trilinos_Belos);
        s.push_back(Solver::Trilinos_ML);
        s.push_back(Solver::Trilinos_Ifpack);
#endif
#if defined(USE_SOLVER_ANI)
        s.push_back(Solver::ANI);
#endif
#if defined(USE_SOLVER_SUPERLU)
        s.push_back(Solver::SUPERLU);
#endif
#if defined(USE_SOLVER_HYPRE)
        s.push_back(Solver::HYPRE_BOOMERAMG);
        s.push_back(Solver::HYPRE_PARASAILS);
        s.push_back(Solver::HYPRE_AMS);
        s.push_back(Solver::HYPRE_PILUT);
        s.push_back(Solver::HYPRE_EUCLID);
#endif
        
#if defined(HAVE_SOLVER_K3BIILU2)
        s.push_back(Solver::K3BIILU2);
#endif

#if defined(HAVE_SOLVER_FCBIILU2)
        s.push_back(Solver::FCBIILU2);
#endif
        return s;
    }

    bool SolverMaster::isSolverAvailable(std::string name)
	{
        try
		{
            SolverInterface *s = SolverMaster::getSolver(name);
            delete s;
            return true;
        }
		catch (...)
		{
            return false;
        }
    }

}

#endif //USE_SOLVER
