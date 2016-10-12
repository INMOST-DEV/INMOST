#include "inmost_solver.h"
#include "inmost_solver_interface.h"

namespace INMOST {

    typedef std::map<std::string, SolverParameter>::const_iterator parameters_map_iterator_t;

    SolverParameters::SolverParameters() {}

    void SolverParameters::SetParameter(std::string name, std::string value) {
        parameters[name] = SolverParameter(value);
    }

    const SolverParameter SolverParameters::GetParameter(std::string name) const {
        parameters_map_iterator_t value = parameters.find(name);
        if (value != parameters.cend()) {
            return value->second;
        } else {
            throw INMOST::SolverUnknownParameter;
        }
    }

}
