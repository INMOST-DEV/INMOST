#include "inmost_solver.h"
#include "inmost_solver_interface.h"

namespace INMOST {

    typedef std::map<std::string, INMOST_DATA_REAL_TYPE>::const_iterator parameters_map_real_iterator_t;
    typedef std::map<std::string, INMOST_DATA_ENUM_TYPE>::const_iterator parameters_map_enum_iterator_t;

    SolverParameters::SolverParameters() {
        this->SetDefaults();
    }

    void SolverParameters::SetDefaults() {
        this->SetParameterEnum("additive_schwartz_overlap", DEFAULT_ADDITIVE_SCHWARTZ_OVERLAP);
        this->SetParameterEnum("maximum_iterations", DEFAULT_MAXIMUM_ITERATIONS);
        this->SetParameterEnum("reorder_nonzero", DEFAULT_PRECONDITIONER_REORDER_NONZEROS);
        this->SetParameterEnum("rescale_iterations", DEFAULT_PRECONDITIONER_RESCALE_ITERS);
        this->SetParameterEnum("condition_estimation", DEFAULT_PRECONDITIONER_CONDITION_ESTIMATION);
        this->SetParameterEnum("adapt_ddpq_tolerance", DEFAULT_PRECONDITIONER_ADAPT_DDPQ_TOLERANCE);
        this->SetParameterEnum("gmres_substeps", DEFAULT_SOLVER_GMRES_SUBSTEPS);

        this->SetParameterReal("absolute_tolerance", DEFAULT_ABSOLUTE_TOLERANCE);
        this->SetParameterReal("relative_tolerance", DEFAULT_RELATIVE_TOLERANCE);
        this->SetParameterReal("divergence_tolerance", DEFAULT_DIVERGENCE_TOLERANCE);
        this->SetParameterReal("drop_tolerance", DEFAULT_PRECONDITIONER_DROP_TOLERANCE);
        this->SetParameterReal("reuse_tolerance", DEFAULT_PRECONDITIONER_REUSE_TOLERANCE);
        this->SetParameterReal("ddpq_tolerance", DEFAULT_PRECONDITIONER_DDPQ_TOLERANCE);
        this->SetParameterReal("fill_level", DEFAULT_PRECONDITIONER_FILL_LEVEL);

    }

    void SolverParameters::SetParameterReal(std::string name, INMOST_DATA_REAL_TYPE value) {
        reals[name] = value;
    }

    void SolverParameters::SetParameterEnum(std::string name, INMOST_DATA_ENUM_TYPE value) {
        enums[name] = value;
    }

    const INMOST_DATA_REAL_TYPE SolverParameters::GetParameterReal(std::string name) const {
        parameters_map_real_iterator_t value = reals.find(name);
        if (value != reals.cend()) {
            return value->second;
        } else {
            throw INMOST::SolverUnknownParameter;
        }
    }

    const INMOST_DATA_ENUM_TYPE SolverParameters::GetParameterEnum(std::string name) const {
        parameters_map_enum_iterator_t value = enums.find(name);
        if (value != enums.cend()) {
            return value->second;
        } else {
            throw INMOST::SolverUnknownParameter;
        }
    }

}
