#ifndef INMOST_SOLVER_MASTER
#define INMOST_SOLVER_MASTER

#include <../Solver/SolverInterface.h>
#if defined(USE_SOLVER)
namespace INMOST {

    class SolverMaster {
    public:
        static SolverInterface *getSolver(std::string name);

        static std::vector<std::string> getAvailableSolvers();

        static bool isSolverAvailable(std::string name);
    };
}
#endif //USE_SOLVER
#endif //INMOST_SOLVER_MASTER
