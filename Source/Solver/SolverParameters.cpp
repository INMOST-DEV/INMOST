#include "inmost_solver.h"
#include "inmost_solver_interface.h"

namespace INMOST {

    const INMOST_DATA_INTEGER_TYPE SolverParameter::integer() const {
        return static_cast<INMOST_DATA_INTEGER_TYPE>(atoi(value.c_str()));
    }

    const INMOST_DATA_ENUM_TYPE SolverParameter::unsigned_integer() const {
        return static_cast<INMOST_DATA_ENUM_TYPE>(atoi(value.c_str()));
    }

    const INMOST_DATA_REAL_TYPE SolverParameter::real() const {
        return static_cast<INMOST_DATA_REAL_TYPE>(atof(value.c_str()));
    }

    const std::string SolverParameter::str() const {
        return std::string(value);
    }

    typedef std::map<std::string, SolverParameter>::const_iterator parameters_map_iterator_t;

    SolverParameters::SolverParameters() {

    }

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
