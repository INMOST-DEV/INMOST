//
// Created by Dmitri Bagaev on 2019-02-20.
//

#ifndef INMOST_TTSP_ANNEALING_H
#define INMOST_TTSP_ANNEALING_H

#include <inmost_ttsp.h>
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

        std::size_t count;
        double      current_value;

        mutable AnnealingUniformDistribution random;

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

        double GetCurrentTemp() const;

        double GetDecrement() const;

        bool IsAllowOscillation() const;

        double GetOscillationTemp() const;

        bool IsStrictBound() const;

        bool IsUseClosest() const;

        double GetNextValue() const;

        double GetRandom() const;

        void SetValue(double value);

        double GetCurrentValue() const;
    };

    class AnnealingOptimizer : public OptimizerInterface {
    private:
        std::size_t                            current_handler_index;
        std::vector<AnnealingParameterHandler> handlers;
        std::vector<double>                    values;
    protected:
        void UpdateSpaceWithLatestResults() override;
    public:
        AnnealingOptimizer(const OptimizationParameters &space, const OptimizerProperties &properties, std::size_t buffer_capacity);

        OptimizationParametersSuggestion Suggest(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                      const OptimizationParameterPoints &,
                                                                                                      void *)> &invoke, void *data) const override;

        virtual ~AnnealingOptimizer();
    };

}


#endif //INMOST_TTSP_ANNEALING_H
