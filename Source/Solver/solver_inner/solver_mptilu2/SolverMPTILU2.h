#ifndef INMOST_SOLVERMPTILU2_H
#define INMOST_SOLVERMPTILU2_H


#include <inmost.h>
#include "solver_mtilu2.hpp"
#include "../SolverInner.h"

namespace INMOST {

    class SolverMPTILU2 : public SolverInner {
    public:
        SolverMPTILU2();

        SolverMPTILU2(const SolverInterface *other);

        virtual void SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner);

        virtual const std::string SolverName() const;

        virtual ~SolverMPTILU2();
    };

}


#endif //INMOST_SOLVERMPTILU2_H
