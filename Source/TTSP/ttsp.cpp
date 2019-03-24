//
// Created by bvdmitri on 04.02.19.
//

#include <inmost_ttsp.h>

#include "optimizers/noop/ttsp_noop.h"
#include "optimizers/bruteforce/ttsp_bruteforce.h"
#include "optimizers/alternating/ttsp_alternating.h"
#include "optimizers/annealing/ttsp_annealing.h"

#if defined(USE_TTSP_LIMBO)

#include "optimizers/bayesian/ttsp_bayesian.h"

#endif

namespace TTSP {

    void OptimizationParameter::swap(OptimizationParameter &left, OptimizationParameter &right) {
        std::swap(left.name, right.name);
        std::swap(left.values, right.values);
        std::swap(left.default_value, right.default_value);
        std::swap(left.type, right.type);
    }

    OptimizationParameter::OptimizationParameter(const std::string &name, const OptimizationParameterRange &range, double step,
                                                 double default_value, OptimizationParameterType type) :
            name(name), values(static_cast<unsigned long>((range.second - range.first) / step)), default_value(default_value), type(type) {
        int    index = 0;
        double value = range.first;
        while (value < range.second) {
            values[index++] = value;
            value += step;
        }
    }

    OptimizationParameter::OptimizationParameter(const std::string &name, const std::vector<double> &values, double default_value, OptimizationParameterType type) :
            name(name), values(values), default_value(default_value), type(type) {}

    OptimizationParameter::OptimizationParameter(const OptimizationParameter &other) :
            name(other.name), values(other.values), default_value(other.default_value), type(other.type) {}

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

    std::size_t OptimizationParameter::GetClosestIndexTo(double to) const noexcept {
        return std::min_element(values.begin(), values.end(), [to](double x, double y) {
            return std::abs(x - to) < std::abs(y - to);
        }) - values.begin();
    }

    std::size_t OptimizationParameter::GetValuesCount() const noexcept {
        return values.size();
    }

    double OptimizationParameter::GetDefaultValue() const noexcept {
        return default_value;
    }

    OptimizationParameterType OptimizationParameter::GetType() const noexcept {
        return type;
    }

    double OptimizationParameter::ExtractValueFromPoint(const OptimizationParameterPoint &point) const noexcept {
        switch (type) {
            case OptimizationParameterType::PARAMETER_TYPE_DEFAULT:
                return point.GetValue();
            case OptimizationParameterType::PARAMETER_TYPE_EXPONENT:
                return std::log10(point.GetValue());
            default:
                return point.GetValue();
        }
    }

    OptimizationParameterPoint::OptimizationParameterPoint(const std::string &name, double value, OptimizationParameterType type) :
            name(name), value(OptimizationParameterPoint::convert(value, type)) {}

    void OptimizationParameterPoint::swap(OptimizationParameterPoint &left, OptimizationParameterPoint &right) {
        std::swap(left.name, right.name);
        std::swap(left.value, right.value);
    }

    double OptimizationParameterPoint::convert(double value, OptimizationParameterType type) {
        switch (type) {
            case OptimizationParameterType::PARAMETER_TYPE_DEFAULT:
                return value;
            case OptimizationParameterType::PARAMETER_TYPE_EXPONENT:
                return std::pow(10, value);
            default:
                return value;
        }
    }

    OptimizationParameterPoint::OptimizationParameterPoint(const OptimizationParameterPoint &other) : name(other.name), value(other.value) {}

    OptimizationParameterPoint::OptimizationParameterPoint(OptimizationParameterPoint &&other) noexcept {
        OptimizationParameterPoint::swap(*this, other);
    }

    OptimizationParameterPoint &OptimizationParameterPoint::operator=(const OptimizationParameterPoint &other) {
        OptimizationParameterPoint tmp(other);
        OptimizationParameterPoint::swap(*this, tmp);
        return *this;
    }

    const std::string &OptimizationParameterPoint::GetName() const noexcept {
        return name;
    }

    double OptimizationParameterPoint::GetValue() const noexcept {
        return value;
    }

    OptimizationParameterPoint::~OptimizationParameterPoint() {}

    OptimizationParametersSuggestion::OptimizationParametersSuggestion(std::size_t changed_index, double changed_value,
                                                                       const OptimizationParameterPoints &before, double metrics_before,
                                                                       const OptimizationParameterPoints &after) :
            changed_index(changed_index), changed_value(changed_value), before(before), metrics_before(metrics_before), after(after) {}

    std::size_t OptimizationParametersSuggestion::GetChangedParameterIndex() const noexcept {
        return changed_index;
    }

    double OptimizationParametersSuggestion::GetChangedValue() const noexcept {
        return changed_value;
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

    std::size_t OptimizationParameters::Size() const {
        return entries.size();
    }

    const OptimizationParameterPoints OptimizationParameters::GetPoints() const noexcept {
        OptimizationParameterPoints points;
        points.reserve(entries.size());

        std::for_each(entries.cbegin(), entries.cend(), [&points](const OptimizationParametersEntry &p) {
            points.emplace_back(OptimizationParameterPoint(p.first.GetName(), p.second, p.first.GetType()));
        });

        return points;
    }

    double OptimizationParameters::GetMetrics() const noexcept {
        return metrics;
    }

    const OptimizationParameterPoints OptimizationParameters::GetPointsWithChangedParameter(const OptimizationParameter &parameter,
                                                                                            double value) const noexcept {
        OptimizationParameterPoints points;
        points.reserve(entries.size());

        std::for_each(entries.cbegin(), entries.cend(), [&points, &parameter, value](const OptimizationParametersEntry &p) {
            points.emplace_back(p.first.GetName() == parameter.GetName() ?
                                OptimizationParameterPoint(p.first.GetName(), value, p.first.GetType()) :
                                OptimizationParameterPoint(p.first.GetName(), p.second, p.first.GetType()));
        });

        return points;
    }

    void OptimizationParameters::Update(const std::vector<double> &update, double metrics) {
        for (int i = 0; i < update.size(); ++i) {
            this->entries[i].second = update[i];
        }

        this->metrics = metrics;
    }

    void OptimizationParameters::Update(std::size_t index, double value, double metrics) {
        this->entries[index].second = value;
        this->metrics               = metrics;
    }

    void OptimizationParameterResult::swap(OptimizationParameterResult &left, OptimizationParameterResult &right) {
        std::swap(left.changed_index, right.changed_index);
        std::swap(left.changed_value, right.changed_value);
        std::swap(left.before, right.before);
        std::swap(left.metrics_before, right.metrics_before);
        std::swap(left.after, right.after);
        std::swap(left.metrics_after, right.metrics_after);
        std::swap(left.is_good, right.is_good);
    }

    OptimizationParameterResult::OptimizationParameterResult(std::size_t changed_index, double changed_value,
                                                             const OptimizationParameterPoints &before, const OptimizationParameterPoints &after,
                                                             double metrics_before, double metrics_after, bool is_good) :
            changed_index(changed_index), changed_value(changed_value),
            before(before), after(after), metrics_before(metrics_before), metrics_after(metrics_after), is_good(is_good) {}

    OptimizationParameterResult::OptimizationParameterResult(const OptimizationParameterResult &other) :
            changed_index(other.changed_index), changed_value(other.changed_value), before(other.before), after(other.after),
            metrics_before(other.metrics_before), metrics_after(other.metrics_after), is_good(other.is_good) {}

    OptimizationParameterResult::OptimizationParameterResult(OptimizationParameterResult &&other) noexcept {
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

    std::size_t OptimizationParameterResult::GetChangedParameterIndex() const noexcept {
        return changed_index;
    }

    /// Getter for new changed value
    double OptimizationParameterResult::GetChangedValue() const noexcept {
        return changed_value;
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

    std::deque<OptimizationParameterResult>::const_reverse_iterator OptimizationParameterResultsBuffer::cbegin() const noexcept {
        return buffer.crbegin();
    }

    const OptimizationParameterResult &OptimizationParameterResultsBuffer::at(std::size_t index) const {
        return *(cbegin() + index);
    }

    std::deque<OptimizationParameterResult>::const_reverse_iterator OptimizationParameterResultsBuffer::cend() const noexcept {
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
        return IsEmpty() ? false : (*cbegin()).IsGood();
    }

    bool OptimizationParameterResultsBuffer::IsSuccessfulResultExist() const noexcept {
        return std::find_if(cbegin(), cend(), [](const OptimizationParameterResult &r) {
            return r.IsGood();
        }) != cend();
    }

    const OptimizationParameterResult &OptimizationParameterResultsBuffer::GetLastSuccessfulResult() const {
        return *std::find_if(cbegin(), cend(), [](const OptimizationParameterResult &r) {
            return r.IsGood();
        });
    }

    const std::vector<OptimizationParameterResult> OptimizationParameterResultsBuffer::GetLastUniqueEntries(std::size_t max) const {
        std::vector<OptimizationParameterResult> unique;

        std::for_each(cbegin(), cend(), [&unique, max](const OptimizationParameterResult &r) {

            if (unique.size() <= max) {

                auto same = std::find_if(unique.cbegin(), unique.cend(), [&r](const OptimizationParameterResult &u) {
                    for (std::size_t i = 0; i < r.GetPointsAfter().size(); ++i) {
                        auto rpi = r.GetPointsAfter().at(i);
                        auto upi = u.GetPointsAfter().at(i);

                        if (std::abs(rpi.GetValue() - upi.GetValue()) > 1e-8) {
                            return false;
                        }
                    }
                    return true;
                });

                if (same == unique.cend()) {
                    unique.push_back(r);
                }

            }

        });

        return unique;
    }

    void OptimizerInterface::RestartWithBestStrategy() {
        auto min = std::min_element(results.cbegin(), results.cbegin() + max_fails, [](const OptimizationParameterResult &l, const OptimizationParameterResult &r) {
            return l.GetMetricsAfter() < r.GetMetricsAfter();
        });

        results.push(OptimizationParameterResult(*min));
        parameters.Update(min->GetChangedParameterIndex(), min->GetChangedValue(), min->GetMetricsAfter());

        this->fails_count = 0;
    }

    bool OptimizerInterface::UpdateSpaceWithLatestResults() {
        if (!results.IsEmpty()) {
            const OptimizationParameterResult &last = results.at(0);
            parameters.Update(last.GetChangedParameterIndex(), last.GetChangedValue(), last.GetMetricsAfter());
        }
        return true;
    }

    OptimizationParametersSuggestion OptimizerInterface::Suggest(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                                      const OptimizationParameterPoints &,
                                                                                                                      void *)> &invoke, void *data) {
        OptimizationAlgorithmSuggestion algorithm = this->AlgorithmMakeSuggestion(invoke, data);

        OptimizationParametersSuggestion suggestion = OptimizationParametersSuggestion(
                algorithm.first,
                algorithm.second,
                parameters.GetPoints(),
                parameters.GetMetrics(),
                parameters.GetPointsWithChangedParameter(parameters.GetParameter(algorithm.first), algorithm.second)
        );

        if (!suggestions.empty()) {
            suggestions.pop_front();
        }

        suggestions.push_front(suggestion);

        return suggestion;
    }

    void OptimizerInterface::SaveResult(std::size_t changed_index, double changed_value,
                                        const OptimizationParameterPoints &before, double metrics_before,
                                        const OptimizationParameterPoints &after, double metrics_after, bool is_good) {
        results.push(OptimizationParameterResult(changed_index, changed_value, before, after, metrics_before, metrics_after, is_good));

        bool is_updated = is_good && this->UpdateSpaceWithLatestResults();
        if (is_updated) {
            this->fails_count = 0;
        } else {
            this->fails_count += 1;
        }

        if (this->fails_count > this->max_fails) {
            switch (restart_strategy) {
                case OptimizerRestartStrategy::RESTART_STRATEGY_NO_RESTART:
                    break;
                case OptimizerRestartStrategy::RESTART_STRATEGY_WITH_BEST:
                    RestartWithBestStrategy();
                    break;
                default:
                    break;
            }
        }
    }

    void OptimizerInterface::SaveResult(std::size_t changed_index, double changed_value,
                                        const OptimizationParameterPoints &after, double metrics_after, bool is_good) {
        SaveResult(changed_index, changed_value, parameters.GetPoints(), parameters.GetMetrics(), after, metrics_after, is_good);
    }

    void OptimizerInterface::SaveResult(const OptimizationParametersSuggestion &suggestion, double metrics, bool is_good) {
        SaveResult(suggestion.GetChangedParameterIndex(), suggestion.GetChangedValue(), suggestion.GetPointsAfter(), metrics, is_good);
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

    const OptimizationParametersSuggestion &OptimizerInterface::GetLastSuggestion() const {
        return suggestions.at(0);
    }

    void OptimizerInterface::SetRestartStrategy(OptimizerRestartStrategy strategy, std::size_t max_fails) noexcept {
        this->restart_strategy = strategy;
        this->max_fails        = max_fails;
    }

    bool Optimizers::IsOptimizerAvailable(const std::string &type) {
        std::vector<std::string> available = Optimizers::GetAvailableOptimizers();
        return std::find(available.cbegin(), available.cend(), type) != available.cend();
    }

    SavedOptimizersMap Optimizers::optimizers = SavedOptimizersMap();

    std::vector<std::string> Optimizers::GetAvailableOptimizers() {
        std::vector<std::string> available;

        available.emplace_back("noop");
        available.emplace_back("bruteforce");
        available.emplace_back("alternating");
        available.emplace_back("annealing");
#if defined(USE_TTSP_LIMBO)
        available.emplace_back("bayesian");
#endif

        return available;
    }

    OptimizerInterface *Optimizers::GetOptimizer(const std::string &type,
                                                 const OptimizationParameters &space, const OptimizerProperties &properties,
                                                 std::size_t buffer_capacity) {
        if (type == "noop") return new NoopOptimizer(space, properties, buffer_capacity);
        if (type == "bruteforce") return new BruteforceOptimizer(space, properties, buffer_capacity);
        if (type == "alternating") return new AlternatingOptimizer(space, properties, buffer_capacity);
        if (type == "annealing") return new AnnealingOptimizer(space, properties, buffer_capacity);
#if defined(USE_TTSP_LIMBO)
        if (type == "bayesian") return new BayesianOptimizer(space, properties, buffer_capacity);
#endif
        return nullptr;
    }

    void Optimizers::SaveOptimizerOrReplace(const std::string &name, const std::string &type, const TTSP::OptimizationParameters &space,
                                            const TTSP::OptimizerProperties &properties, std::size_t buffer_capacity) {

        OptimizerInterface *created = Optimizers::GetOptimizer(type, space, properties, buffer_capacity);

        auto optimizer = Optimizers::optimizers.find(name);
        if (optimizer != Optimizers::optimizers.end()) {
            delete (*optimizer).second;
        }
        Optimizers::optimizers[name] = created;
    }

    OptimizerInterface *Optimizers::GetSavedOptimizer(const std::string &name) {
        auto optimizer = Optimizers::optimizers.find(name);
        if (optimizer == Optimizers::optimizers.end()) {
            return nullptr;
        } else {
            return (*optimizer).second;
        }
    }
}
