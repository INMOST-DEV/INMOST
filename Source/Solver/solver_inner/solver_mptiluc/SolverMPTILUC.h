#ifndef INMOST_SOLVERMPTILUC_H
#define INMOST_SOLVERMPTILUC_H

#include <inmost.h>
#include "solver_mtiluc2.hpp"
#include "../SolverInner.h"

namespace INMOST {

    class SolverMPTILUC : public SolverInner {
    public:
        SolverMPTILUC();

        SolverMPTILUC(const SolverInterface *other);

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner);

        virtual bool Solve(Sparse::Vector &RHS, Sparse::Vector &SOL);

        virtual const std::string SolverName() const;

        virtual ~SolverMPTILUC();
    };

}


#endif //INMOST_SOLVERMPTILUC_H
