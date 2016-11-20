#ifndef INMOST_SOLVERDDPQILUC2_H
#define INMOST_SOLVERDDPQILUC2_H


#include <inmost.h>
#include "solver_ddpqiluc2.hpp"
#include "../SolverInner.h"

namespace INMOST {

    class SolverDDPQILUC2 : public SolverInner {
    public:
        SolverDDPQILUC2();

        SolverDDPQILUC2(const SolverInterface *other);

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner);

        virtual const std::string SolverName() const;

        virtual ~SolverDDPQILUC2();
    };

}


#endif //INMOST_SOLVERDDPQILUC2_H
