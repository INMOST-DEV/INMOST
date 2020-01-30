#ifndef INMOST_SOLVERSUPERLU_H
#define INMOST_SOLVERSUPERLU_H

#include "../Solver/SolverInterface.h"

#if defined(USE_SOLVER_SUPERLU_DIST)
#include "superlu_ddefs.h"
#else
#include "superlu/slu_ddefs.h"
#endif

namespace INMOST {

    class SolverSUPERLU : public SolverInterface {
    private:
#if defined(USE_SOLVER_SUPERLU_DIST)
		LUstruct_t LUstruct;
		SOLVEstruct_t SOLVEstruct;
		ScalePermstruct_t ScalePermstruct;
		gridinfo_t grid;
		superlu_dist_options_t options_;
		SuperLUStat_t stat_;
		int g_size;
#else //USE_SOLVER_SUPERLU_DIST
        int * perm_r;
        int * perm_c;
        int * remap;
        superlu_options_t options;
        SuperLUStat_t stat;
#endif //USE_SOLVER_SUPERLU_DIST
		int a_size, info;
        SuperMatrix A, L, U;
    public:
        SolverSUPERLU();

        virtual SolverInterface *Copy(const SolverInterface *other);

        virtual void Assign(const SolverInterface *other);

        virtual void Setup(int *argc, char ***argv, SolverParameters &p);

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner);

        virtual bool Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL);

        virtual bool Clear();

        virtual bool isMatrixSet();

        virtual std::string GetParameter(std::string name) const;

        virtual void SetParameter(std::string name, std::string value);

        virtual INMOST_DATA_ENUM_TYPE Iterations() const;

        virtual INMOST_DATA_REAL_TYPE Residual() const;

        virtual const std::string ReturnReason() const;

        virtual const std::string SolverName() const;

        virtual void Finalize();

        virtual ~SolverSUPERLU();

    };

}


#endif //INMOST_SOLVERSUPERLU_H
