//
// Created by Dmitri Bagaev on 2019-02-20.
//

#include "ttsp_annealing.h"

namespace TTSP {

    template<typename T>
    int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }


    AnnealingUniformDistribution::AnnealingUniformDistribution() : distribution(0.0, 1.0) {
        int rank, size;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        unsigned int seed = static_cast<unsigned int>(time(NULL));

        if (rank == 0) {
            for (int i = 1; i < size; i++) {
                MPI_Send(&seed, 1, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD);
            }
        } else {
            MPI_Status status;
            MPI_Recv(&seed, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &status);
        }

        generator.seed(seed);
    }

    double AnnealingUniformDistribution::next() {
        return distribution(generator);
    }

    double                               AnnealingParameterHandler::DEFAULT_TEMP0             = 0.0010037;
    double                               AnnealingParameterHandler::DEFAULT_DECREMENT         = 2;
    bool                                 AnnealingParameterHandler::DEFAULT_ALLOW_OSCILLATION = false;
    double                               AnnealingParameterHandler::DEFAULT_OSCILLATION_TEMP  = 9.76563e-10;
    bool                                 AnnealingParameterHandler::DEFAULT_STRICT_BOUND      = true;
    bool                                 AnnealingParameterHandler::DEFAULT_USE_CLOSEST       = false;

    AnnealingParameterHandler::AnnealingParameterHandler(const OptimizationParameter &parameter, const OptimizerInterface &optimizer) :
            count(0), parameter(parameter), value(parameter.GetDefaultValue()),
            temp0(AnnealingParameterHandler::DEFAULT_TEMP0),
            decrement(AnnealingParameterHandler::DEFAULT_DECREMENT),
            allow_oscillation(AnnealingParameterHandler::DEFAULT_ALLOW_OSCILLATION),
            oscillation_temp(AnnealingParameterHandler::DEFAULT_OSCILLATION_TEMP),
            strict_bound(AnnealingParameterHandler::DEFAULT_STRICT_BOUND),
            use_closest(AnnealingParameterHandler::DEFAULT_USE_CLOSEST) {

        const std::string &temp0_property = parameter.GetName() + ":temp0";
        if (optimizer.HasProperty(temp0_property)) {
            temp0 = std::atof(optimizer.GetProperty(temp0_property).c_str());
        }

        const std::string &decrement_property = parameter.GetName() + ":decrement";
        if (optimizer.HasProperty(decrement_property)) {
            decrement = std::atof(optimizer.GetProperty(decrement_property).c_str());
        }

        const std::string &allow_oscillation_property = parameter.GetName() + ":allow_oscillation";
        if (optimizer.HasProperty(decrement_property)) {
            const std::string &allow = optimizer.GetProperty(decrement_property);
            if (allow == "true") {
                allow_oscillation = true;
            } else if (allow == "false") {
                allow_oscillation = false;
            }
        }

        const std::string &oscillation_temp_property = parameter.GetName() + ":oscillation_temp";
        if (optimizer.HasProperty(oscillation_temp_property)) {
            oscillation_temp = std::atof(optimizer.GetProperty(oscillation_temp_property).c_str());
        }

        const std::string &strict_bound_property = parameter.GetName() + ":strict_bound";
        if (optimizer.HasProperty(strict_bound_property)) {
            const std::string &strict = optimizer.GetProperty(strict_bound_property);
            if (strict == "true") {
                strict_bound = true;
            } else if (strict == "false") {
                strict_bound = false;
            }
        }

        const std::string &use_closest_property = parameter.GetName() + ":use_closest";
        if (optimizer.HasProperty(use_closest_property)) {
            const std::string &closest = optimizer.GetProperty(use_closest_property);
            if (closest == "true") {
                use_closest = true;
            } else if (closest == "false") {
                use_closest = false;
            }
        }
    }

    AnnealingParameterHandler::AnnealingParameterHandler(const AnnealingParameterHandler &other) :
            count(other.count), value(other.value),
            parameter(other.parameter), temp0(other.temp0), decrement(other.decrement),
            allow_oscillation(other.allow_oscillation), oscillation_temp(other.oscillation_temp), strict_bound(other.strict_bound),
            use_closest(other.use_closest) {}

    const OptimizationParameter &AnnealingParameterHandler::GetParameter() const {
        return parameter;
    }

    double AnnealingParameterHandler::GetTemp0() const {
        return temp0;
    }

    double AnnealingParameterHandler::GetCurrentTemp() const {
        return temp0 * std::exp(-decrement * std::pow(count, 1.0 / 2.0));
    }

    double AnnealingParameterHandler::GetDecrement() const {
        return decrement;
    }

    bool AnnealingParameterHandler::IsAllowOscillation() const {
        return allow_oscillation;
    }

    double AnnealingParameterHandler::GetOscillationTemp() const {
        return oscillation_temp;
    }

    bool AnnealingParameterHandler::IsStrictBound() const {
        return strict_bound;
    }

    bool AnnealingParameterHandler::IsUseClosest() const {
        return use_closest;
    }

    double AnnealingParameterHandler::GetNextValue() const {

        double current = GetCurrentValue();
        double temp    = GetCurrentTemp();
        double alpha   = random.next();
        double z       = sgn(alpha - 1.0 / 2.0) * temp * (std::pow(1.0 + 1.0 / temp, std::abs(2 * alpha - 1)) - 1);

        double a     = parameter.GetMinimalValue();
        double b     = parameter.GetMaximumValue();
        double bound = (b - a);

        if (strict_bound) {
            double aa = current - a;
            double bb = b - current;

            if (bb < aa) {
                bound = 2 * bb;
            } else {
                bound = 2 * aa;
            }
        }

        double next = current + z * bound;
        while (!((next >= a) && (next <= b))) {
            alpha = random.next();
            z     = sgn(alpha - 1.0 / 2.0) * temp * (std::pow(1.0 + 1.0 / temp, std::abs(2 * alpha - 1)) - 1);
            next  = current + z * bound;
        }

        if (use_closest) {
            next = parameter.GetClosestTo(next);
        }

        return next;
    }

    double AnnealingParameterHandler::GetRandom() const {
        return random.next();
    }

    double AnnealingParameterHandler::GetCurrentValue() const {
        return value;
    }

    void AnnealingParameterHandler::SetValue(double value) {
        this->count += 1;
        this->value = value;
    }

    AnnealingOptimizer::AnnealingOptimizer(const OptimizationParameters &space, const OptimizerProperties &properties, std::size_t buffer_capacity) :
            OptimizerInterface(space, properties, buffer_capacity), current_handler_index(0) {
        const OptimizationParameterEntries &parameters = space.GetParameterEntries();
        handlers.reserve(parameters.size());
        std::for_each(parameters.cbegin(), parameters.cend(), [this](const OptimizationParametersEntry &entry) {
            handlers.emplace_back(AnnealingParameterHandler(entry.first, *this));
        });
    }

    OptimizationAlgorithmSuggestion AnnealingOptimizer::AlgorithmMakeSuggestion(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                                                     const OptimizationParameterPoints &,
                                                                                                                                     void *)> &invoke, void *data) const {
        return std::make_pair(current_handler_index, handlers.at(current_handler_index).GetNextValue());
    }

    bool AnnealingOptimizer::UpdateSpaceWithLatestResults() {
        AnnealingParameterHandler         &h    = handlers.at(current_handler_index);
        const OptimizationParameterResult &last = results.at(0);

        double temp    = h.GetCurrentTemp();
        double delta_e = last.GetMetricsAfter() - last.GetMetricsBefore();
        double et      = 1.0 / (1.0 + std::exp(delta_e / temp));
        double alpha   = h.GetRandom();

        bool is_updated = false;

        if (last.IsGood() && (last.GetMetricsBefore() < 0.0 || ((delta_e < 0.0) || alpha < et))) {
            double update_value = last.GetChangedValue();
            h.SetValue(update_value);
            parameters.Update(current_handler_index, update_value, last.GetMetricsAfter());
            is_updated = true;
        }

        current_handler_index = (current_handler_index + 1) % (handlers.size());

        return is_updated;
    }

    AnnealingOptimizer::~AnnealingOptimizer() {}
}
