#ifndef INMOST_SOLVERTRILINOSIFPACK_H
#define INMOST_SOLVERTRILINOSIFPACK_H

#include <inmost.h>
#include "../SolverTrilinos.h"

namespace INMOST {

    class SolverTrilinosIfpack: public SolverTrilinos {
    public:
        SolverTrilinosIfpack(SolverParameters &parameters);

        SolverTrilinosIfpack(const SolverInterface *other);

        virtual bool Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL);

        virtual const std::string SolverName() const;

        virtual ~SolverTrilinosIfpack();
    };

}


#endif //INMOST_SOLVERTRILINOSIFPACK_H
