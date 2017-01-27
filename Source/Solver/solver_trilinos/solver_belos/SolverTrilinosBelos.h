#ifndef INMOST_SOLVERTRILINOSBELOS_H
#define INMOST_SOLVERTRILINOSBELOS_H

#include <inmost.h>
#include "../SolverTrilinos.h"

namespace INMOST {

    class SolverTrilinosBelos : public SolverTrilinos {
    public:
        SolverTrilinosBelos();

        SolverTrilinosBelos(const SolverInterface *other);

        virtual bool Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL);

        virtual const std::string SolverName() const;

        virtual ~SolverTrilinosBelos();
    };

}


#endif //INMOST_SOLVERTRILINOSBELOS_H
