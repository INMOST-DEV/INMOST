#ifndef INMOST_SOLVER_MASTER
#define INMOST_SOLVER_MASTER

#include <Source/Solver/SolverInterface.h>

namespace INMOST {

    class SolverMaster {
    public:
        static SolverInterface *getSolver(std::string name);

        static std::vector<std::string> getAvailableSolvers();

        static bool isSolverAvailable(std::string name);
    };
}

#endif //INMOST_SOLVER_MASTER
