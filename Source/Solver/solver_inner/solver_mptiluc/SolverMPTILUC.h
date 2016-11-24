#ifndef INMOST_SOLVERMPTILUC_H
#define INMOST_SOLVERMPTILUC_H

#include <inmost.h>
#include "solver_mtiluc2.hpp"
#include "../SolverInner.h"

namespace INMOST {

    class SolverMPTILUC : public SolverInner {
        INMOST_DATA_ENUM_TYPE rescale_iterations, condition_estimation, schwartz_overlap, gmres_substeps, reorder_nnz, fill_level;
        INMOST_DATA_REAL_TYPE drop_tolerance, reuse_tolerance;
    public:
        SolverMPTILUC();

        SolverMPTILUC(const SolverInterface *other);

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner);

        virtual void SetParameter(std::string name, std::string value);

        virtual const std::string SolverName() const;

        virtual ~SolverMPTILUC();
    };

}


#endif //INMOST_SOLVERMPTILUC_H
