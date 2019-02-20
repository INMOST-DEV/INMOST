//
// Created by Dmitri Bagaev on 2019-02-20.
//

#include "ttsp_annealing.h"

namespace TTSP {

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
            parameter(parameter),
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

    AnnealingOptimizer::AnnealingOptimizer(const OptimizationParametersSpace &space) : OptimizerInterface(space, 10), current_handler_index(0) {
        const OptimizationParameters &parameters = space.GetParameters();
        handlers.reserve(parameters.size());
        std::for_each(parameters.begin(), parameters.end(), [this](const OptimizationParametersEntry &entry) {
            handlers.emplace_back(AnnealingParameterHandler(entry.first, *this));
        });
    }

    OptimizationParameterPoints AnnealingOptimizer::MakeOptimizationIteration(INMOST::Solver &solver, INMOST::Sparse::Matrix &matrix,
                                                                              INMOST::Sparse::Vector &RHS) {
        return TTSP::OptimizationParameterPoints();
    }

    AnnealingOptimizer::~AnnealingOptimizer() {}
}
