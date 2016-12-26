#ifndef INMOST_SOLVERAZTEC_H
#define INMOST_SOLVERAZTEC_H

#include "SolverTrilinosAztec.h"
#include "../SolverTrilinos.h"

namespace INMOST {

    class SolverTrilinosAztec: public SolverTrilinos {
    public:
        SolverTrilinosAztec();

        SolverTrilinosAztec(const SolverInterface *other);

        virtual bool Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL);

        virtual const std::string SolverName() const;

        virtual ~SolverTrilinosAztec();
    };

}

#endif //INMOST_SOLVERAZTEC_H
