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

    OptimizationParametersSuggestion::OptimizationParametersSuggestion(const OptimizationParameter &changed,
                                                                       const OptimizationParameterPoints &before, double metrics_before,
                                                                       const OptimizationParameterPoints &after) :
            changed(changed), before(before), metrics_before(metrics_before), after(after) {}

    const OptimizationParameter &OptimizationParametersSuggestion::GetChangedParameter() const noexcept {
        return changed;
    }

    const OptimizationParameterPoints &OptimizationParametersSuggestion::GetPointsBefore() const noexcept {
        return before;
    }

    double OptimizationParametersSuggestion::GetMetricsBefore() const noexcept {
        return metrics_before;
    }

    const OptimizationParameterPoints &OptimizationParametersSuggestion::GetPointsAfter() const noexcept {
        return after;
    }

    void OptimizationParameters::swap(OptimizationParameters &left, OptimizationParameters &right) {
        std::swap(left.entries, right.entries);
        std::swap(left.metrics, right.metrics);
    }

    OptimizationParameters::OptimizationParameters(const OptimizationParameterEntries &entries, double metrics) : entries(entries), metrics(metrics) {}

    OptimizationParameters::OptimizationParameters(const OptimizationParameters &other) : entries(other.entries), metrics(other.metrics) {}

    OptimizationParameters::OptimizationParameters(OptimizationParameters &&other) noexcept {
        OptimizationParameters::swap(*this, other);
    }

    OptimizationParameters &OptimizationParameters::operator=(const OptimizationParameters &other) {
        OptimizationParameters tmp(other);
        OptimizationParameters::swap(*this, tmp);
        return *this;
    }

    const OptimizationParameterEntries &OptimizationParameters::GetParameterEntries() const noexcept {
        return entries;
    }

    const OptimizationParametersEntry &OptimizationParameters::GetParameterEntry(std::size_t index) const {
        return entries.at(index);
    }

    const OptimizationParameter &OptimizationParameters::GetParameter(std::size_t index) const {
        return GetParameterEntry(index).first;
    }

    const OptimizationParameterPoints OptimizationParameters::GetPoints() const noexcept {
        OptimizationParameterPoints points(entries.size());
        std::transform(entries.begin(), entries.end(), points.begin(), [](const OptimizationParametersEntry &p) {
            return std::make_pair(p.first.GetName(), p.second);
        });
        return points;
    }

    double OptimizationParameters::GetMetrics() const noexcept {
        return metrics;
    }

    const OptimizationParameterPoints OptimizationParameters::GetPointsWithChangedParameter(const OptimizationParameter &parameter,
                                                                                            double value) const noexcept {
        OptimizationParameterPoints points(entries.size());
        std::transform(entries.begin(), entries.end(), points.begin(), [&parameter, value](const OptimizationParametersEntry &p) {
            return p.first.GetName() == parameter.GetName() ? std::make_pair(p.first.GetName(), value) : std::make_pair(p.first.GetName(), p.second);
        });
        return points;
    }

    void OptimizationParameters::Update(const OptimizationParameterPoints &update, double metrics) {
        for (int i = 0; i < update.size(); ++i) {
            this->entries[i].second = update[i].second;
        }

        this->metrics = metrics;
    }

    void OptimizationParameters::Update(std::size_t index, double value, double metrics) {
        this->entries[index].second = value;
        this->metrics               = metrics;
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

    void OptimizerInterface::UpdateSpaceWithLatestResults() {
        if (!results.IsEmpty()) {
            const OptimizationParameterResult &result = (*results.begin());
            parameters.Update(result.GetPointsAfter(), result.GetMetricsAfter());
        }
    }

    void OptimizerInterface::SaveResult(const OptimizationParameter &changed,
                                        const OptimizationParameterPoints &before, double metrics_before,
                                        const OptimizationParameterPoints &after, double metrics_after, bool is_good) {
        results.push(OptimizationParameterResult(changed, before, after, metrics_before, metrics_after, is_good));
        if (is_good) {
            this->UpdateSpaceWithLatestResults();
        }
    }

    void OptimizerInterface::SaveResult(const OptimizationParameter &changed, const OptimizationParameterPoints &after, double metrics_after, bool is_good) {
        SaveResult(changed, parameters.GetPoints(), parameters.GetMetrics(), after, metrics_after, is_good);
    }

    const OptimizationParameterResultsBuffer &OptimizerInterface::GetResults() const noexcept {
        return results;
    };

    const OptimizationParameterPoints OptimizerInterface::GetPoints() const noexcept {
        return parameters.GetPoints();
    }

    void OptimizerInterface::SetVerbosityLevel(OptimizerVerbosityLevel level) noexcept {
        verbosity = level;
    }

    OptimizerVerbosityLevel OptimizerInterface::GetVerbosityLevel() const noexcept {
        return verbosity;
    }

    bool OptimizerInterface::HasProperty(const std::string &name) const noexcept {
        return properties.find(name) != properties.end();
    }

    const std::string &OptimizerInterface::GetProperty(const std::string &name) const {
        return properties.at(name);
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
                                                         const OptimizationParameters &space, const OptimizerProperties &properties,
                                                         std::size_t buffer_capacity) {
        if (type == "noop") return new NoopOptimizer(space, properties, buffer_capacity);
        if (type == "bruteforce") return new BruteforceOptimizer(space, properties, buffer_capacity);
        if (type == "alternating") return new AlternatingOptimizer(space, properties, buffer_capacity);
        if (type == "annealing") return new AnnealingOptimizer(space, properties, buffer_capacity);
        return nullptr;
    }
}
