//
// Created by Dmitri Bagaev on 2019-03-22.
//

#ifndef INMOST_BAYESIAN_H
#define INMOST_BAYESIAN_H


#include <inmost_optimizer.h>
#include <random>

namespace INMOST {

    class BayesianUniformDistribution {
    private:
        std::mt19937                     generator;
        std::uniform_real_distribution<> distribution;
    public:
        explicit BayesianUniformDistribution();

        double next();

        unsigned int operator()(unsigned int n);
    };

    class BayesianOptimizer : public OptimizerInterface {
    private:
        static unsigned int DEFAULT_UNIQUE_POINTS_MAX_COUNT;
        static unsigned int DEFAULT_UNIQUE_POINTS_RANDOM_COUNT;
        static unsigned int DEFAULT_INITIAL_ITERATIONS_COUNT;
        static double       DEFAULT_INITIAL_ITERATIONS_RADIUS;
        static double       DEFAULT_MAX_JUMP_BARRIER;

        unsigned int unique_points_max_count;
        unsigned int unique_points_random_count;
        unsigned int initial_iterations_count;
        double       initial_iterations_radius;
        double       max_jump_barrier;

        mutable BayesianUniformDistribution random;

    protected:
        bool UpdateSpaceWithLatestResults() override;

        SuggestionChangedParameters AlgorithmMakeSuggestion() const override;

    public:
        BayesianOptimizer(const std::string &name, const OptimizationParameters &space, const OptimizerProperties &properties, std::size_t buffer_capacity);

        virtual ~BayesianOptimizer();
    };

}


#endif //INMOST_BAYESIAN_H
