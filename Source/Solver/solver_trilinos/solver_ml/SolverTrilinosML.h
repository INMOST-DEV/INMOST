#ifndef INMOST_SOLVERTRILINOSML_H
#define INMOST_SOLVERTRILINOSML_H

#include <inmost.h>
#include "../SolverTrilinos.h"

namespace INMOST {

    class SolverTrilinosML: public SolverTrilinos {
    public:
        SolverTrilinosML();

        SolverTrilinosML(const SolverInterface *other);

        virtual bool Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL);

        virtual const std::string SolverName() const;

        virtual ~SolverTrilinosML();
    };

}


#endif //INMOST_SOLVERTRILINOSML_H
