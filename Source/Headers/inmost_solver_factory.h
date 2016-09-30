#ifndef INMOST_SOLVER_FACTORY
#define INMOST_SOLVER_FACTORY

namespace INMOST {

    struct SolverBaseFactory {
        virtual SolverInterface *create() = 0;
        virtual SolverInterface *copy(const SolverInterface *other) = 0;
        virtual ~SolverBaseFactory() {};
    };

    template<class C>
    struct SolverCreateFactory : SolverBaseFactory {
        SolverInterface *create() {
            return new C;
        };

        SolverInterface *copy(const SolverInterface *other) {
            return new C(other);
        };
    };

    class SolverFactory {
    private:
        static std::map<std::string, SolverBaseFactory *> solvers;
    public:
        template<class T>
        static void registerSolver(std::string name) {
            solvers.insert(std::make_pair(name, new SolverCreateFactory<T>));
        };
        static SolverInterface *getSolver(std::string name);
        static SolverInterface *copySolver(const SolverInterface *other);
        static std::vector<std::string> getAvailableSolvers();
        static bool isSolverAvailable(std::string name);
    };

}

#endif //INMOST_SOLVER_FACTORY
