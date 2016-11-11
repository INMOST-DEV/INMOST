#include "inmost_solver_parameters.h"
#include "inmost_solver.h"

namespace INMOST {

    SolverParameters::SolverParameters() {}

    SolverParameters::SolverParameters(std::string solverName, std::string prefix) : solverName(solverName), prefix(prefix) {}

    SolverParameters::SolverParameters(const SolverParameters &other) : prefix(other.prefix), parameters(other.parameters) {}

    void SolverParameters::set(std::string name, std::string value) {
        parameters[name] = value;
    }

    void SolverParameters::require(std::string name, std::string value) {
        if (parameters.find(name) == parameters.end()) {
            parameters[name] = value;
        }
    }

    std::string SolverParameters::getSolverName() const {
        return solverName;
    }

    std::string SolverParameters::getPrefix() const {
        return prefix;
    }

    typedef std::vector<SolverParameters>::iterator parameters_iterator_t;
    std::vector<SolverParameters> GlobalSolversParameters::parameters = std::vector<SolverParameters>();

    SolverParameters &GlobalSolversParameters::getSolverParameters(std::string solverName, std::string prefix) {
        for (parameters_iterator_t i = parameters.begin(); i < parameters.end(); i++) {
            SolverParameters& p = *i;
            if (p.getSolverName() == solverName && p.getPrefix() == prefix) return p;
        }
        SolverParameters new_p = SolverParameters(solverName, prefix);
        parameters.push_back(new_p);
        return (parameters[parameters.size() - 1]);
    }

    void GlobalSolversParameters::Initialize(std::string databasePath) {
        std::vector<std::string> solvers = Solver::getAvailableSolvers();
        for (std::vector<std::string>::iterator name = solvers.begin(); name < solvers.end(); name++) {
            //global[*name] = std::map<std::string, SolverParameters>();
        }
    }
}
