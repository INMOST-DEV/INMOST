//
// Created by Dmitri Bagaev on 2019-03-22.
//

#ifndef INMOST_TTSP_BAYESIAN_H
#define INMOST_TTSP_BAYESIAN_H


#include <inmost_ttsp.h>
#include <random>

namespace TTSP {

    class BayesianUniformDistribution {
    private:
        std::mt19937                     generator;
        std::uniform_real_distribution<> distribution;
    public:
        explicit BayesianUniformDistribution();

        double next();
    };

    class BayesianOptimizer : public OptimizerInterface {
    private:
        static unsigned int DEFAULT_INITIAL_ITERATIONS_COUNT;
        static double       DEFAULT_INITIAL_ITERATIONS_RADIUS;

        unsigned int initial_iterations_count;
        double       initial_iterations_radius;

        mutable BayesianUniformDistribution random;

    protected:
        OptimizationAlgorithmSuggestion AlgorithmMakeSuggestion(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                                     const OptimizationParameterPoints &,
                                                                                                                     void *)> &invoke, void *data) const override;

    public:
        BayesianOptimizer(const OptimizationParameters &space, const OptimizerProperties &properties, std::size_t buffer_capacity);

        virtual ~BayesianOptimizer();
    };

}


#endif //INMOST_TTSP_BAYESIAN_H
