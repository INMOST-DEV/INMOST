#ifndef INMOST_SOLVERMLMPTILUC_H
#define INMOST_SOLVERMLMPTILUC_H

#include <inmost.h>
#include "solver_mlmtiluc2.hpp"
#include "../SolverInner.h"

namespace INMOST
{

    class SolverMLMPTILUC : public SolverInner
	{
        INMOST_DATA_ENUM_TYPE rescale_iterations, condition_estimation, schwartz_overlap, gmres_substeps, reorder_nnz, fill_level, verbosity;
        INMOST_DATA_REAL_TYPE drop_tolerance, reuse_tolerance, pivot_condition, pivot_diag;
    public:
        SolverMLMPTILUC();

        SolverMLMPTILUC(const SolverInterface *other);

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner);

        virtual void SetParameter(std::string name, std::string value);

        virtual const std::string SolverName() const;

        virtual ~SolverMLMPTILUC();
    };

}


#endif //INMOST_SOLVERMLMPTILUC_H
