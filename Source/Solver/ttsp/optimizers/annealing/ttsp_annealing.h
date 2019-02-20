//
// Created by Dmitri Bagaev on 2019-02-20.
//

#ifndef INMOST_TTSP_ANNEALING_H
#define INMOST_TTSP_ANNEALING_H

#include <Source/Solver/ttsp/ttsp.h>
#include <random>

namespace TTSP {

    class AnnealingUniformDistribution {
    private:
        std::mt19937                     generator;
        std::uniform_real_distribution<> distribution;
    public:
        explicit AnnealingUniformDistribution();

        double next();
    };

    enum AnnealingParameterOptimizationType {
        STRICT,
        EXPONENT
    };

    class AnnealingParameterHandler {
    private:
        static AnnealingParameterOptimizationType DEFAULT_OPTIMIZATION_TYPE;
        static double                             DEFAULT_TEMP0;
        static double                             DEFAULT_DECREMENT;
        static bool                               DEFAULT_ALLOW_OSCILLATION;
        static double                             DEFAULT_OSCILLATION_TEMP;
        static bool                               DEFAULT_STRICT_BOUND;
        static bool                               DEFAULT_USE_CLOSEST;

        const OptimizationParameter        &parameter;
        AnnealingParameterOptimizationType optimization_type;
        double                             temp0;
        double                             decrement;
        bool                               allow_oscillation;
        double                             oscillation_temp;
        bool                               strict_bound;
        bool                               use_closest;
    public:
        explicit AnnealingParameterHandler(const OptimizationParameter &parameter, const OptimizerInterface &optimizer);

        AnnealingParameterHandler(const AnnealingParameterHandler &other);

        const OptimizationParameter &GetParameter() const;

        AnnealingParameterOptimizationType GetOptimizationType() const;

        double GetTemp0() const;

        double GetDecrement() const;

        bool IsAllowOscillation() const;

        double GetOscillationTemp() const;

        bool IsStrictBound() const;

        bool IsUseClosest() const;
    };

    class AnnealingOptimizer : public OptimizerInterface {
    private:
        AnnealingUniformDistribution           random;
        std::size_t                            current_handler_index;
        std::vector<AnnealingParameterHandler> handlers;
    public:
        AnnealingOptimizer(const OptimizationParametersSpace &space);

        OptimizationParameterPoints MakeOptimizationIteration(INMOST::Solver &solver, INMOST::Sparse::Matrix &matrix,
                                                              INMOST::Sparse::Vector &RHS) override;

        virtual ~AnnealingOptimizer();
    };

}


#endif //INMOST_TTSP_ANNEALING_H
