#ifndef INMOST_SOLVERILU2_H
#define INMOST_SOLVERILU2_H

#include <inmost.h>
#include "solver_ilu2.hpp"
#include "../SolverInner.h"

namespace INMOST {

    class SolverILU2 : public SolverInner {
        INMOST_DATA_ENUM_TYPE rescale_iterations, schwartz_overlap, gmres_substeps, reorder_nnz, fill_level, verbosity;
        INMOST_DATA_REAL_TYPE drop_tolerance, reuse_tolerance;
    public:
        SolverILU2();

        SolverILU2(const SolverInterface *other);

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner);

        virtual void SetParameter(std::string name, std::string value);

        virtual const std::string SolverName() const;

        virtual ~SolverILU2();
    };

}

#endif
