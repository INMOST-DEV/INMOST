//
// Created by bvdmitri on 04.02.19.
//

#include <Source/Solver/ttsp/optimizers/noop/ttsp_noop.h>
#include <Source/Solver/ttsp/optimizers/bruteforce/ttsp_bruteforce.h>
#include <Source/Solver/ttsp/optimizers/alternating/ttsp_alternating.h>
#include <Source/Solver/ttsp/optimizers/annealing/ttsp_annealing.h>
#include "ttsp.h"

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
        int    index = 0;
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

    const std::string &OptimizationParameter::GetName() const noexcept {
        return name;
    }

    const std::vector<double> &OptimizationParameter::GetValues() const noexcept {
        return values;
    }

    double OptimizationParameter::GetMinimalValue() const noexcept {
        return values.at(0);
    }

    double OptimizationParameter::GetMaximumValue() const noexcept {
        return values.at(values.size() - 1);
    }

    double OptimizationParameter::GetClosestTo(double to) const noexcept {
        return *std::min_element(values.begin(), values.end(), [to](double x, double y) {
            return std::abs(x - to) < std::abs(y - to);
        });
    }

    std::size_t OptimizationParameter::GetValuesCount() const noexcept {
        return values.size();
    }

    double OptimizationParameter::GetDefaultValue() const noexcept {
        return default_value;
    }

    OptimizationParametersSuggestion::OptimizationParametersSuggestion(const OptimizationParameter &changed, const OptimizationParameterPoints &before,
                                                                       const OptimizationParameterPoints &after) :
            changed(changed), before(before), after(after) {}

    const OptimizationParameter &OptimizationParametersSuggestion::GetChangedParameter() const noexcept {
        return changed;
    }

    const OptimizationParameterPoints &OptimizationParametersSuggestion::GetPointsBefore() const noexcept {
        return before;
    }

    const OptimizationParameterPoints &OptimizationParametersSuggestion::GetPointsAfter() const noexcept {
        return after;
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

    bool OptimizationParametersSpace::isSolverNameMatch(const std::string &solver_name) const noexcept {
        return this->solver_name == solver_name;
    }

    bool OptimizationParametersSpace::isSolverPrefixMatch(const std::string &solver_prefix) const noexcept {
        return this->solver_prefix == solver_prefix;
    }

    const std::string &OptimizationParametersSpace::GetSolverName() const noexcept {
        return solver_name;
    }

    const std::string &OptimizationParametersSpace::GetSolverPrefix() const noexcept {
        return solver_prefix;
    }

    const OptimizationParameters &OptimizationParametersSpace::GetParameters() const noexcept {
        return parameters;
    }

    const OptimizationParameterPoints OptimizationParametersSpace::GetPoints() const noexcept {
        OptimizationParameterPoints points(parameters.size());
        std::transform(parameters.begin(), parameters.end(), points.begin(), [](const OptimizationParametersEntry &p) {
            return std::make_pair(p.first.GetName(), p.second);
        });
        return points;
    }

    const OptimizationParameterPoints OptimizationParametersSpace::GetPointsWithChangedParameter(const OptimizationParameter &parameter,
                                                                                                 double value) const noexcept {
        OptimizationParameterPoints points(parameters.size());
        std::transform(parameters.begin(), parameters.end(), points.begin(), [&parameter, value](const OptimizationParametersEntry &p) {
            return p.first.GetName() == parameter.GetName() ? std::make_pair(p.first.GetName(), value) : std::make_pair(p.first.GetName(), p.second);
        });
        return points;
    }

    void OptimizationParametersSpace::Update(const OptimizationParameterPoints &update) {
        for (int i = 0; i < update.size(); ++i) {
            parameters[i].second = update[i].second;
        }
    }

    void OptimizationParameterResult::swap(OptimizationParameterResult &left, OptimizationParameterResult &right) {
        std::swap(left.before, right.before);
        std::swap(left.metrics_before, right.metrics_before);
        std::swap(left.after, right.after);
        std::swap(left.metrics_after, right.metrics_after);
        std::swap(left.is_good, right.is_good);
    }

    OptimizationParameterResult::OptimizationParameterResult(const OptimizationParameter &changed, const OptimizationParameterPoints &before,
                                                             const OptimizationParameterPoints &after,
                                                             double metrics_before, double metrics_after, bool is_good) :
            changed(changed), before(before), after(after), metrics_before(metrics_before), metrics_after(metrics_after), is_good(is_good) {}

    OptimizationParameterResult::OptimizationParameterResult(const OptimizationParameterResult &other) :
            changed(other.changed), before(other.before), after(other.after),
            metrics_before(other.metrics_before), metrics_after(other.metrics_after), is_good(other.is_good) {}

    OptimizationParameterResult::OptimizationParameterResult(OptimizationParameterResult &&other) noexcept : changed(std::move(other).changed) {
        OptimizationParameterResult::swap(*this, other);
    }

    OptimizationParameterResult &OptimizationParameterResult::operator=(const OptimizationParameterResult &other) {
        OptimizationParameterResult tmp(other);
        OptimizationParameterResult::swap(*this, tmp);
        return *this;
    }

    const OptimizationParameterPoints &OptimizationParameterResult::GetPointsBefore() const noexcept {
        return before;
    }

    const OptimizationParameterPoints &OptimizationParameterResult::GetPointsAfter() const noexcept {
        return after;
    }

    double OptimizationParameterResult::GetMetricsBefore() const noexcept {
        return metrics_before;
    }

    double OptimizationParameterResult::GetMetricsAfter() const noexcept {
        return metrics_after;
    }

    bool OptimizationParameterResult::IsGood() const noexcept {
        return is_good;
    }


    void OptimizationParameterResultsBuffer::swap(OptimizationParameterResultsBuffer &left,
                                                  OptimizationParameterResultsBuffer &right) {
        std::swap(left.buffer, right.buffer);
        std::swap(left.capacity, right.capacity);
    }

    OptimizationParameterResultsBuffer::OptimizationParameterResultsBuffer(std::size_t capacity) : capacity(capacity) {}

    OptimizationParameterResultsBuffer::OptimizationParameterResultsBuffer(const OptimizationParameterResultsBuffer &other) :
            buffer(other.buffer), capacity(other.capacity) {}

    OptimizationParameterResultsBuffer::OptimizationParameterResultsBuffer(OptimizationParameterResultsBuffer &&other) noexcept {
        OptimizationParameterResultsBuffer::swap(*this, other);
    }

    OptimizationParameterResultsBuffer &OptimizationParameterResultsBuffer::operator=(const OptimizationParameterResultsBuffer &other) {
        OptimizationParameterResultsBuffer tmp(other);
        OptimizationParameterResultsBuffer::swap(*this, tmp);
        return *this;
    }

    std::deque<OptimizationParameterResult>::const_reverse_iterator OptimizationParameterResultsBuffer::begin() const noexcept {
        return buffer.crbegin();
    }

    const OptimizationParameterResult &OptimizationParameterResultsBuffer::at(std::size_t index) const {
        return *(begin() + index);
    }

    std::deque<OptimizationParameterResult>::const_reverse_iterator OptimizationParameterResultsBuffer::end() const noexcept {
        return buffer.crend();
    }

    void OptimizationParameterResultsBuffer::push(const OptimizationParameterResult &result) {
        if (buffer.size() == capacity) {
            buffer.pop_front();
        }
        buffer.push_back(result);
    }

    std::size_t OptimizationParameterResultsBuffer::size() const noexcept {
        return buffer.size();
    }

    bool OptimizationParameterResultsBuffer::IsEmpty() const noexcept {
        return buffer.size() == 0;
    }

    bool OptimizationParameterResultsBuffer::IsLastResultSuccessful() const noexcept {
        return IsEmpty() ? false : (*begin()).IsGood();
    }

    bool OptimizationParameterResultsBuffer::IsSuccessfulResultExist() const noexcept {
        return std::find_if(begin(), end(), [](const OptimizationParameterResult &r) {
            return r.IsGood();
        }) != end();
    }

    const OptimizationParameterResult &OptimizationParameterResultsBuffer::GetLastSuccessfulResult() const {
        return *std::find_if(begin(), end(), [](const OptimizationParameterResult &r) {
            return r.IsGood();
        });
    }

    void OptimizerInterface::SaveResult(const OptimizationParameter &changed,
                                        const OptimizationParameterPoints &before, double metrics_before,
                                        const OptimizationParameterPoints &after, double metrics_after, bool is_good) {
        results.push(OptimizationParameterResult(changed, before, after, metrics_before, metrics_after, is_good));
    }

    void OptimizerInterface::SaveResult(const OptimizationParameter &changed,
                                        const OptimizationParameterPoints &before,
                                        const OptimizationParameterPoints &after,
                                        double metrics_after, bool is_good) {
        results.push(
                OptimizationParameterResult(changed, before, after, results.IsEmpty() ? -1.0 : (*results.begin()).GetMetricsAfter(), metrics_after, is_good)
        );
    }

    const OptimizationParameterResultsBuffer &OptimizerInterface::GetResults() const noexcept {
        return results;
    };

    void OptimizerInterface::SetVerbosityLevel(OptimizerVerbosityLevel level) noexcept {
        verbosity = level;
    }

    OptimizerVerbosityLevel OptimizerInterface::GetVerbosityLevel() const noexcept {
        return verbosity;
    }

    void OptimizerInterface::SetProperty(const std::string &name, const std::string &value) {
        properties[name] = value;
    }

    bool OptimizerInterface::HasProperty(const std::string &name) const noexcept {
        return properties.find(name) != properties.end();
    }

    const std::string &OptimizerInterface::GetProperty(const std::string &name) const {
        return properties.at(name);
    }

    void OptimizerInterface::UpdateSpacePoints(const OptimizationParameterPoints &update) {
        space.Update(update);
    }

    bool OptimizerInterface::IsOptimizerAvailable(const std::string &type) {
        std::vector<std::string> available = OptimizerInterface::GetAvailableOptimizers();
        return std::find(available.begin(), available.end(), type) != available.end();
    }

    std::vector<std::string> OptimizerInterface::GetAvailableOptimizers() {
        std::vector<std::string> available;

        available.emplace_back("noop");
        available.emplace_back("bruteforce");
        available.emplace_back("alternating");
        available.emplace_back("annealing");

        return available;
    }

    OptimizerInterface *OptimizerInterface::GetOptimizer(const std::string &type,
                                                         const OptimizationParametersSpace &space, const OptimizerProperties &properties,
                                                         std::size_t buffer_capacity) {
        if (type == "noop") return new NoopOptimizer(space, properties, buffer_capacity);
        if (type == "bruteforce") return new BruteforceOptimizer(space, properties, buffer_capacity);
        if (type == "alternating") return new AlternatingOptimizer(space, properties, buffer_capacity);
        if (type == "annealing") return new AnnealingOptimizer(space, properties, buffer_capacity);
        return nullptr;
    }
}
