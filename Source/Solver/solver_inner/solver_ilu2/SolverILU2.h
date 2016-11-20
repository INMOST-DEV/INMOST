#ifndef INMOST_SOLVERILU2_H
#define INMOST_SOLVERILU2_H

#include <inmost.h>
#include "solver_ilu2.hpp"
#include "../SolverInner.h"

namespace INMOST {

    class SolverILU2 : public SolverInner {
    public:
        SolverILU2();

        SolverILU2(const SolverInterface *other);

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner);

        virtual const std::string SolverName() const;

        virtual ~SolverILU2();
    };

}

#endif