#ifndef INMOST_SOLVERDDPQILUC2_H
#define INMOST_SOLVERDDPQILUC2_H


#include <inmost.h>
#include "solver_ddpqiluc2.hpp"
#include "../SolverInner.h"

namespace INMOST {

    class SolverDDPQILUC2 : public SolverInner {
        INMOST_DATA_ENUM_TYPE rescale_iterations, condition_estimation, adapt_ddpq_tolerance, schwartz_overlap, gmres_substeps, reorder_nnz;
        INMOST_DATA_REAL_TYPE ddpq_tolerance, drop_tolerance, reuse_tolerance;
    public:
    public:
        SolverDDPQILUC2();

        SolverDDPQILUC2(const SolverInterface *other);

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner);

        virtual void SetParameter(std::string name, std::string value);

        virtual const std::string SolverName() const;

        virtual ~SolverDDPQILUC2();
    };

}


#endif //INMOST_SOLVERDDPQILUC2_H
