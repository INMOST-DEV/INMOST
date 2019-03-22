//
// Created by Dmitri Bagaev on 2019-03-22.
//

#ifndef INMOST_TTSP_BAYESIAN_H
#define INMOST_TTSP_BAYESIAN_H


#include <inmost_ttsp.h>

namespace TTSP {

    class BayesianOptimizer : public OptimizerInterface {
    private:
        mutable std::size_t current_iteration;
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
