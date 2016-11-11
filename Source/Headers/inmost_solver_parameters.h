#ifndef INMOST_SOLVER_PARAMETERS_H
#define INMOST_SOLVER_PARAMETERS_H

#include <string>
#include "inmost_utils.h"
#include "inmost_common.h"

namespace INMOST {

#if defined(USE_SOLVER)

    class SolverParameters {
    private:
        std::map<std::string, std::string> parameters;

        std::string prefix;
        std::string solverName;
    public:
        SolverParameters();

        SolverParameters(std::string solverName, std::string prefix);

        SolverParameters(const SolverParameters &other);

        std::string getSolverName() const;
        std::string getPrefix() const;

        void require(std::string name, std::string value);
        void set(std::string name, std::string value);

        template<typename T>
        T get(std::string name) {
            return from_string<T>(parameters[name]);
        }

    };

    class GlobalSolversParameters {
    private:
        static std::vector<SolverParameters> parameters;

        GlobalSolversParameters() {};

        static void Initialize(std::string databasePath);
    public:
        static SolverParameters& getSolverParameters(std::string solverName, std::string prefix);
    };

#endif

}

#endif //INMOST_SOLVER_PARAMETERS_H
