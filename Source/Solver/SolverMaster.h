#ifndef INMOST_SOLVER_MASTER
#define INMOST_SOLVER_MASTER

#include <inmost_solver.h>
#include <inmost_solver_interface.h>

namespace INMOST {

    struct SolverBaseFactory {
        virtual SolverInterface *create(SolverParameters &parameters) = 0;

        virtual SolverInterface *copy(const SolverInterface *other) = 0;

        virtual ~SolverBaseFactory() {};
    };

    template<class C>
    struct SolverCreateFactory : SolverBaseFactory {
        SolverInterface *create(SolverParameters &parameters) {
            return new C(parameters);
        };

        SolverInterface *copy(const SolverInterface *other) {
            return new C(other);
        };
    };

    class SolverMaster {
    private:
        static std::map<std::string, SolverBaseFactory *> solvers;

        template<class T>
        static void registerSolver(std::string name) {
            solvers.insert(std::make_pair(name, new SolverCreateFactory<T>));
        };

        static SolverInterface *getSolver(std::string name, SolverParameters &parameters);

        static SolverInterface *copySolver(const SolverInterface *other);

        static std::vector<std::string> getAvailableSolvers();

        static bool isSolverAvailable(std::string name);

        friend Solver::Solver(std::string solverName, std::string prefix, INMOST_MPI_Comm _comm);

        friend Solver::Solver(const Solver &other);

        friend void Solver::Initialize(int *argc, char ***argv, const char *database);

        friend bool Solver::isSolverAvailable(std::string name);

        friend std::vector<std::string> Solver::getAvailableSolvers();
    };

    typedef std::map<std::string, SolverBaseFactory *>::iterator solvers_map_iterator_t;
}

#endif //INMOST_SOLVER_MASTER
