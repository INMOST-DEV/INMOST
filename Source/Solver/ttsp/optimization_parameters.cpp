//
// Created by bvdmitri on 31.01.19.
//

#include "optimization_parameters.h"

#include <algorithm>

namespace TTSP {

    void OptimizationParameter::swap(OptimizationParameter &left, OptimizationParameter &right) {
        std::swap(left.name, right.name);
        std::swap(left.values, right.values);
        std::swap(left.default_value, right.default_value);
    }

    OptimizationParameter::OptimizationParameter(const std::string &name,
                                                 const OptimizationParameterRange &range,
                                                 double step, double default_value) :
            name(name),
            values(static_cast<unsigned long>((range.second - range.first) / step)),
            default_value(default_value) {
        int index = 0;
        double value = range.first;
        while (value < range.second) {
            values[index++] = value;
            value += step;
        }
    }

    OptimizationParameter::OptimizationParameter(const std::string &name, const std::vector<double> &values,
                                                 double default_value) :
            name(name), values(values), default_value(default_value) {}

    OptimizationParameter::OptimizationParameter(const OptimizationParameter &other) :
            name(other.name), values(other.values), default_value(other.default_value) {}

    OptimizationParameter::OptimizationParameter(OptimizationParameter &&other) noexcept {
        OptimizationParameter::swap(*this, other);
    }

    OptimizationParameter &TTSP::OptimizationParameter::operator=(const OptimizationParameter &other) {
        OptimizationParameter tmp(other);
        OptimizationParameter::swap(*this, tmp);
        return *this;

    }

    const std::string &OptimizationParameter::GetName() const {
        return name;
    }

    const std::vector<double> &OptimizationParameter::GetValues() const {
        return values;
    }

    const double &OptimizationParameter::GetDefaultValue() const {
        return default_value;
    }

    void OptimizationParametersSpace::swap(OptimizationParametersSpace &left, OptimizationParametersSpace &right) {
        std::swap(left.solver_name, right.solver_name);
        std::swap(left.solver_prefix, right.solver_prefix);
        std::swap(left.parameters, right.parameters);
    }

    OptimizationParametersSpace::OptimizationParametersSpace(const std::string &solver_name,
                                                             const std::string &solver_prefix,
                                                             const OptimizationParameters &parameters) :
            solver_name(solver_name), solver_prefix(solver_prefix), parameters(parameters) {}

    OptimizationParametersSpace::OptimizationParametersSpace(const OptimizationParametersSpace &other) :
            solver_name(other.solver_name), solver_prefix(other.solver_prefix), parameters(other.parameters) {}

    OptimizationParametersSpace::OptimizationParametersSpace(OptimizationParametersSpace &&other) noexcept {
        OptimizationParametersSpace::swap(*this, other);
    }

    OptimizationParametersSpace &OptimizationParametersSpace::operator=(const OptimizationParametersSpace &other) {
        OptimizationParametersSpace tmp(other);
        OptimizationParametersSpace::swap(*this, tmp);
        return *this;
    }

    bool OptimizationParametersSpace::isSolverNameMatch(const std::string &solver_name) const {
        return this->solver_name == solver_name;
    }

    bool OptimizationParametersSpace::isSolverPrefixMatch(const std::string &solver_prefix) const {
        return this->solver_prefix == solver_prefix;
    }

    const std::string &OptimizationParametersSpace::GetSolverName() const {
        return solver_name;
    }

    const std::string &OptimizationParametersSpace::GetSolverPrefix() const {
        return solver_prefix;
    }

    const OptimizationParameters &OptimizationParametersSpace::GetParameters() const {
        return parameters;
    }

    const OptimizationParameterPoints OptimizationParametersSpace::GetPoints() const {
        OptimizationParameterPoints points(parameters.size());
        std::transform(parameters.begin(), parameters.end(), points.begin(), [](const OptimizationParametersEntry &p) {
            return std::make_pair(p.first.GetName(), p.second);
        });
        return points;
    }

    void OptimizationParametersSpace::Update(const OptimizationParameterPoints &update) {
        for (int i = 0; i < update.size(); ++i) {
            parameters[i].second = update[i].second;
        }
    }
}


