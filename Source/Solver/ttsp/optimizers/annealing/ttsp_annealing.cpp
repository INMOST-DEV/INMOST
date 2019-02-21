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

    AnnealingParameterOptimizationType   AnnealingParameterHandler::DEFAULT_OPTIMIZATION_TYPE = STRICT;
    double                               AnnealingParameterHandler::DEFAULT_TEMP0             = 0.0010037;
    double                               AnnealingParameterHandler::DEFAULT_DECREMENT         = 2;
    bool                                 AnnealingParameterHandler::DEFAULT_ALLOW_OSCILLATION = false;
    double                               AnnealingParameterHandler::DEFAULT_OSCILLATION_TEMP  = 9.76563e-10;
    bool                                 AnnealingParameterHandler::DEFAULT_STRICT_BOUND      = true;
    bool                                 AnnealingParameterHandler::DEFAULT_USE_CLOSEST       = true;

    AnnealingParameterHandler::AnnealingParameterHandler(const OptimizationParameter &parameter, const OptimizerInterface &optimizer) :
            count(0), parameter(parameter), current_value(parameter.GetDefaultValue()),
            optimization_type(AnnealingParameterHandler::DEFAULT_OPTIMIZATION_TYPE),
            temp0(AnnealingParameterHandler::DEFAULT_TEMP0),
            decrement(AnnealingParameterHandler::DEFAULT_DECREMENT),
            allow_oscillation(AnnealingParameterHandler::DEFAULT_ALLOW_OSCILLATION),
            oscillation_temp(AnnealingParameterHandler::DEFAULT_OSCILLATION_TEMP),
            strict_bound(AnnealingParameterHandler::DEFAULT_STRICT_BOUND),
            use_closest(AnnealingParameterHandler::DEFAULT_USE_CLOSEST) {

        const std::string &optimization_type_property = parameter.GetName() + ":optimization_type";
        if (optimizer.HasProperty(optimization_type_property)) {
            const std::string &type = optimizer.GetProperty(optimization_type_property);
            if (type == "strict") {
                optimization_type = STRICT;
            } else if (type == "exponent") {
                optimization_type = EXPONENT;
            }
        }

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
            count(other.count), current_value(other.current_value),
            parameter(other.parameter), optimization_type(other.optimization_type), temp0(other.temp0), decrement(other.decrement),
            allow_oscillation(other.allow_oscillation), oscillation_temp(other.oscillation_temp), strict_bound(other.strict_bound),
            use_closest(other.use_closest) {}

    const OptimizationParameter &AnnealingParameterHandler::GetParameter() const {
        return parameter;
    }

    AnnealingParameterOptimizationType AnnealingParameterHandler::GetOptimizationType() const {
        return optimization_type;
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

        if (optimization_type == EXPONENT) {
            return std::pow(10, -next);
        } else {
            return next;
        }
    }

    double AnnealingParameterHandler::GetRandom() const {
        return random.next();
    }

    void AnnealingParameterHandler::SetValue(double value) {
        current_value = value;
        count += 1;
    }

    double AnnealingParameterHandler::GetCurrentValue() const {
        return current_value;
    }

    AnnealingOptimizer::AnnealingOptimizer(const OptimizationParametersSpace &space, const OptimizerProperties &properties, std::size_t buffer_capacity) :
            OptimizerInterface(space, properties, buffer_capacity), current_handler_index(0), values(space.GetParameters().size()) {
        const OptimizationParameters &parameters = space.GetParameters();
        handlers.reserve(parameters.size());
        std::for_each(parameters.begin(), parameters.end(), [this](const OptimizationParametersEntry &entry) {
            handlers.emplace_back(AnnealingParameterHandler(entry.first, *this));
        });
        std::transform(parameters.begin(), parameters.end(), values.begin(), [](const OptimizationParametersEntry &entry) {
            return entry.first.GetDefaultValue();
        });
    }

    OptimizationParameterPoints AnnealingOptimizer::MakeOptimizationIteration(INMOST::Solver &solver, INMOST::Sparse::Matrix &matrix,
                                                                              INMOST::Sparse::Vector &RHS) {

        OptimizationParameterPoints points(space.GetParameters().size());

        if (results.size() < 2) {
            if (results.size() == 0) {
                std::transform(handlers.begin(), handlers.end(), points.begin(), [](const AnnealingParameterHandler &h) {
                    return std::make_pair(h.GetParameter().GetName(), h.GetCurrentValue());
                });
            } else if (results.size() == 1) {
                int i = 0;
                std::transform(handlers.begin(), handlers.end(), points.begin(), [this, &i](const AnnealingParameterHandler &h) {
                    return std::make_pair(h.GetParameter().GetName(), i++ == current_handler_index ? h.GetNextValue() : h.GetCurrentValue());
                });
            }
        } else {
            AnnealingParameterHandler         &h           = handlers.at(current_handler_index);
            const OptimizationParameterResult &last        = results.at(0);
            const OptimizationParameterResult &before_last = results.at(1);

            double temp    = h.GetCurrentTemp();
            double delta_e = last.GetTime() - before_last.GetTime();
            double et      = 1.0 / (1.0 + std::exp(delta_e / temp));
            //double h = std::exp(-delta_e / temp);
            double alpha   = h.GetRandom();


            if (last.IsSolved() && ((delta_e < 0.0) || alpha < et)) {
                double update_value = last.GetPoints().at(current_handler_index).second;
                h.SetValue(update_value);
            }

            current_handler_index = (current_handler_index + 1) % (handlers.size());

            int i = 0;
            std::transform(handlers.begin(), handlers.end(), points.begin(), [this, &i](AnnealingParameterHandler &h) {
                return std::make_pair(h.GetParameter().GetName(), i++ == current_handler_index ? h.GetNextValue() : h.GetCurrentValue());
            });
        }

        return points;
    }

    AnnealingOptimizer::~AnnealingOptimizer() {}
}
