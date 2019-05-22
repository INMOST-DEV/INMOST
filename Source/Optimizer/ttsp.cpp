//
// Created by bvdmitri on 04.02.19.
//

#include <inmost_optimizer.h>

#include "Source/Optimizer/optimizers/noop/noop.h"
#include "Source/Optimizer/optimizers/bruteforce/bruteforce.h"
#include "Source/Optimizer/optimizers/alternating/alternating.h"
#include "Source/Optimizer/optimizers/annealing/annealing.h"

#if defined(USE_OPTIMIZER_BAYESIAN)

#include "Source/Optimizer/optimizers/bayesian/bayesian.h"

#endif

namespace INMOST {

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

    OptimizationParameter &INMOST::OptimizationParameter::operator=(const OptimizationParameter &other) {
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

    double OptimizationParameter::GetMinStrictBound(double v, double modifier) const noexcept {
        double min = GetMinimalValue();
        double max = GetMaximumValue();

        double aa = v - min;
        double bb = max - v;

        double bound = 0.0;

        if (bb < aa) {
            bound = 2 * bb;
        } else {
            bound = 2 * aa;
        }

        return std::min(v - bound * modifier, max - (max - min) * modifier);
    }

    double OptimizationParameter::GetMaxStrictBound(double v, double modifier) const noexcept {
        double min = GetMinimalValue();
        double max = GetMaximumValue();

        double aa = v - min;
        double bb = max - v;

        double bound = 0.0;

        if (bb < aa) {
            bound = 2 * bb;
        } else {
            bound = 2 * aa;
        }

        return std::max(v + bound * modifier, min + (max - min) * modifier);
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

    void SuggestionChangedParameter::swap(SuggestionChangedParameter &left, SuggestionChangedParameter &right) {
        std::swap(left.index, right.index);
        std::swap(left.value, right.value);
    }

    SuggestionChangedParameter::SuggestionChangedParameter(std::size_t index, double value) : index(index), value(value) {}

    SuggestionChangedParameter::SuggestionChangedParameter(const SuggestionChangedParameter &other) : index(other.index), value(other.value) {}

    SuggestionChangedParameter::SuggestionChangedParameter(SuggestionChangedParameter &&other) {
        SuggestionChangedParameter::swap(*this, other);
    }

    SuggestionChangedParameter &SuggestionChangedParameter::operator=(const SuggestionChangedParameter &other) {
        SuggestionChangedParameter tmp(other);
        SuggestionChangedParameter::swap(*this, tmp);
        return *this;
    }

    std::size_t SuggestionChangedParameter::GetIndex() const noexcept {
        return index;
    }

    double SuggestionChangedParameter::GetValue() const noexcept {
        return value;
    }

    SuggestionChangedParameter::~SuggestionChangedParameter() {}

    OptimizationParametersSuggestion::OptimizationParametersSuggestion(const SuggestionChangedParameters &changed,
                                                                       const OptimizationParameterPoints &before, double metrics_before,
                                                                       const OptimizationParameterPoints &after) :
            changed(changed), before(before), metrics_before(metrics_before), after(after) {}

    const SuggestionChangedParameters &OptimizationParametersSuggestion::GetChangedParameters() const noexcept {
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

    const OptimizationParameterPoints OptimizationParameters::GetPointsWithChangedParameters(const SuggestionChangedParameters &changed) const noexcept {
        OptimizationParameterPoints points;
        points.reserve(entries.size());

        std::size_t index = 0;
        std::for_each(entries.cbegin(), entries.cend(), [&index, &points, &changed](const OptimizationParametersEntry &e) {

            auto found = (std::find_if(changed.cbegin(), changed.cend(), [index](const SuggestionChangedParameter &c) {
                return index == c.GetIndex();
            }));

            points.emplace_back(found == changed.cend() ?
                                OptimizationParameterPoint(e.first.GetName(), e.second, e.first.GetType()) :
                                OptimizationParameterPoint(e.first.GetName(), (*found).GetValue(), e.first.GetType())
            );

            index += 1;

        });

        return points;
    }

    void OptimizationParameters::Update(const SuggestionChangedParameters &changed, double metrics) {
        std::for_each(changed.cbegin(), changed.cend(), [this](const SuggestionChangedParameter &c) {
            this->entries[c.GetIndex()].second = c.GetValue();
        });

        this->metrics = metrics;
    }

    void OptimizationParameters::Update(std::size_t index, double value, double metrics) {
        this->entries[index].second = value;
        this->metrics               = metrics;
    }

    void OptimizationParameterResult::swap(OptimizationParameterResult &left, OptimizationParameterResult &right) {
        std::swap(left.changed, right.changed);
        std::swap(left.before, right.before);
        std::swap(left.metrics_before, right.metrics_before);
        std::swap(left.after, right.after);
        std::swap(left.metrics_after, right.metrics_after);
        std::swap(left.is_good, right.is_good);
    }

    OptimizationParameterResult::OptimizationParameterResult(const SuggestionChangedParameters &changed,
                                                             const OptimizationParameterPoints &before, const OptimizationParameterPoints &after,
                                                             double metrics_before, double metrics_after, bool is_good) :
            changed(changed), before(before), after(after), metrics_before(metrics_before), metrics_after(metrics_after), is_good(is_good) {}

    OptimizationParameterResult::OptimizationParameterResult(const OptimizationParameterResult &other) :
            changed(other.changed), before(other.before), after(other.after),
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

    const SuggestionChangedParameters &OptimizationParameterResult::GetChangedParameters() const noexcept {
        return changed;
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
        parameters.Update(min->GetChangedParameters(), min->GetMetricsAfter());

        this->fails_count = 0;
    }

    bool OptimizerInterface::UpdateSpaceWithLatestResults() {
        if (!results.IsEmpty()) {
            const OptimizationParameterResult &last = results.at(0);
            parameters.Update(last.GetChangedParameters(), last.GetMetricsAfter());
        }
        return true;
    }

    OptimizationParametersSuggestion OptimizerInterface::Suggest(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                                      const OptimizationParameterPoints &,
                                                                                                                      void *)> &invoke, void *data) {
        SuggestionChangedParameters changed = this->AlgorithmMakeSuggestion(invoke, data);

        OptimizationParametersSuggestion suggestion = OptimizationParametersSuggestion(
                changed,
                parameters.GetPoints(),
                parameters.GetMetrics(),
                parameters.GetPointsWithChangedParameters(changed)
        );

        if (!suggestions.empty()) {
            suggestions.pop_front();
        }

        suggestions.push_front(suggestion);

        if (verbosity >= OptimizerVerbosityLevel::Level2 && mpi_rank == 0) {
            std::cout << "[Optimizer:" << name << ":Suggestion]";

            std::for_each(changed.cbegin(), changed.cend(), [this](const SuggestionChangedParameter &c) {
                std::cout << "\t" << parameters.GetParameter(c.GetIndex()).GetName() << "=" << c.GetValue();
            });

            std::cout << std::endl;

        }

        return suggestion;
    }

    void OptimizerInterface::SaveResult(const SuggestionChangedParameters &changed,
                                        const OptimizationParameterPoints &before, double metrics_before,
                                        const OptimizationParameterPoints &after, double metrics_after, bool is_good) {
        results.push(OptimizationParameterResult(changed, before, after, metrics_before, metrics_after, is_good));

        if (verbosity >= OptimizerVerbosityLevel::Level1 && mpi_rank == 0) {
            std::cout << "[Optimizer:" << name << ":Result]";

            std::size_t i = 0;

            for (const OptimizationParametersEntry &entry : parameters.GetParameterEntries()) {
                std::cout << "\t" << "before." << entry.first.GetName() << "=" << before.at(i++).GetValue();
            }
            std::cout << "\tbefore.metrics=" << metrics_before;

            i = 0;
            for (const OptimizationParametersEntry &entry : parameters.GetParameterEntries()) {
                std::cout << "\t" << "after." << entry.first.GetName() << "=" << after.at(i++).GetValue();
            }
            std::cout << "\tafter.metrics=" << metrics_after << std::endl;
        }

        if (verbosity >= INMOST::OptimizerVerbosityLevel::Level3 && mpi_rank == 0) {
            std::cout << "[Optimizer:" << name << ":ResultsBuffer]" << std::endl;

            int index = 1;
            std::for_each(results.cbegin(), results.cend(), [&index](const INMOST::OptimizationParameterResult &result) {
                std::cout << "\t" << index++ << "\t" << " [";

                const INMOST::OptimizationParameterPoints &pbefore = result.GetPointsBefore();
                std::for_each(pbefore.begin(), pbefore.end(), [](const INMOST::OptimizationParameterPoint &point) {
                    std::cout << " " << point.GetName() << "=" << point.GetValue() << " ";
                });

                std::cout << "] -> [";

                const INMOST::OptimizationParameterPoints &pafter = result.GetPointsAfter();
                std::for_each(pafter.begin(), pafter.end(), [](const INMOST::OptimizationParameterPoint &point) {
                    std::cout << " " << point.GetName() << "=" << point.GetValue() << " ";
                });

                std::cout << "]\t(" << result.GetMetricsBefore() << " -> " << result.GetMetricsAfter() << ")" << std::endl;
            });
        }

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

    void OptimizerInterface::SaveResult(const SuggestionChangedParameters &changed,
                                        const OptimizationParameterPoints &after, double metrics_after, bool is_good) {
        SaveResult(changed, parameters.GetPoints(), parameters.GetMetrics(), after, metrics_after, is_good);
    }

    void OptimizerInterface::SaveResult(const OptimizationParametersSuggestion &suggestion, double metrics, bool is_good) {
        SaveResult(suggestion.GetChangedParameters(), suggestion.GetPointsAfter(), metrics, is_good);
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

    const std::string &OptimizerInterface::GetName() const noexcept {
        return name;
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
#if defined(USE_OPTIMIZER_BAYESIAN)
        available.emplace_back("bayesian");
#endif

        return available;
    }

    OptimizerInterface *Optimizers::GetOptimizer(const std::string &name, const std::string &type,
                                                 const OptimizationParameters &parameters, const OptimizerProperties &properties,
                                                 std::size_t buffer_capacity) {
        if (type == "noop") return new NoopOptimizer(name, parameters, properties, buffer_capacity);
        if (type == "bruteforce") return new BruteforceOptimizer(name, parameters, properties, buffer_capacity);
        if (type == "alternating") return new AlternatingOptimizer(name, parameters, properties, buffer_capacity);
        if (type == "annealing") return new AnnealingOptimizer(name, parameters, properties, buffer_capacity);
#if defined(USE_OPTIMIZER_BAYESIAN)
        if (type == "bayesian") return new BayesianOptimizer(name, parameters, properties, buffer_capacity);
#endif
        return nullptr;
    }

    void Optimizers::SaveOptimizerOrReplace(const std::string &name, const std::string &type, const INMOST::OptimizationParameters &parameters,
                                            const INMOST::OptimizerProperties &properties, std::size_t buffer_capacity) {

        if (parameters.Size() == 0) {
            std::cout << "Cannot create optimizer with empty parameters" << std::endl;
            return;
        }

        OptimizerInterface *created = Optimizers::GetOptimizer(name, type, parameters, properties, buffer_capacity);

        auto optimizer = Optimizers::optimizers.find(name);
        if (optimizer != Optimizers::optimizers.end()) {
            delete (*optimizer).second;
        }
        Optimizers::optimizers[name] = created;
    }

    void Optimizers::SaveOptimizerOrReplace(const std::string &name, const std::string &type, const OptimizerProperties &properties, std::size_t buffer_capacity) {
        SaveOptimizerOrReplace(name, type, OptimizersConfiguration::GetGlobalOptimizationParameters(), properties, buffer_capacity);
    }

    OptimizerInterface *Optimizers::GetSavedOptimizer(const std::string &name) {
        auto optimizer = Optimizers::optimizers.find(name);
        if (optimizer == Optimizers::optimizers.end()) {
            return nullptr;
        } else {
            return (*optimizer).second;
        }
    }

    std::vector<OptimizationParameter> OptimizersConfiguration::global_parameters;

    OptimizationParameter OptimizersConfiguration::ParseOptimizationParameter(const std::string &line) {
        std::istringstream line_stream(line);

        std::string type;
        std::string name;
        double      default_value;

        line_stream >> type >> name >> default_value;

        std::vector<double> values;

        while (!line_stream.eof()) {
            double tmp;
            line_stream >> tmp;
            values.push_back(tmp);
        }

        bool is_raw_values = true;

        if (type.find("range", 0) != std::string::npos) {
            is_raw_values = false;
            std::cout << "range " << std::endl;
        }

        OptimizationParameterType ptype = OptimizationParameterType::PARAMETER_TYPE_DEFAULT;

        if (type.find("exponent", 0) != std::string::npos) {
            ptype = OptimizationParameterType::PARAMETER_TYPE_EXPONENT;
        }

        if (is_raw_values) {
            return OptimizationParameter(name, values, default_value, ptype);
        } else {
            return OptimizationParameter(name, std::make_pair(values.at(0), values.at(1)), values.at(2), default_value, ptype);
        }
    }

    void OptimizersConfiguration::FromFile(const std::string &path) {
        std::ifstream file;
        file.open(path);
        if (!file.bad()) {
            FromStream(file);
        } else {
            std::cout << "[OptimizersConfiguration] Cannot open file: " << path << std::endl;
        }
    }

    void OptimizersConfiguration::FromStream(std::istream &stream) {
        OptimizerConfigurationReaderState state = OptimizerConfigurationReaderState::READ_GLOBAL_PARAMETER;

        std::string line;
        while (std::getline(stream, line)) {
            // TODO more features
            switch (state) {
                case OptimizerConfigurationReaderState::READ_GLOBAL_PARAMETER:
                    OptimizersConfiguration::global_parameters.emplace_back(OptimizersConfiguration::ParseOptimizationParameter(line));
                    break;
                default:
                    std::cout << "[OptimizersConfiguration] Unreachable state in optimizers configuration reader" << std::endl;
            }
        }
    }

    OptimizationParameters OptimizersConfiguration::GetGlobalOptimizationParameters() {
        std::vector<OptimizationParametersEntry> entries;
        std::for_each(OptimizersConfiguration::global_parameters.cbegin(), OptimizersConfiguration::global_parameters.cend(), [&entries](const OptimizationParameter &p) {
            entries.emplace_back(std::make_pair(OptimizationParameter(p), p.GetDefaultValue()));
        });
        return OptimizationParameters(entries, -1);
    }
}
