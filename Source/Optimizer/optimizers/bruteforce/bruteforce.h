//
// Created by bvdmitri on 31.01.19.
//

#ifndef INMOST_BRUTEFORCE_H
#define INMOST_BRUTEFORCE_H

#include <inmost_optimizer.h>

namespace INMOST {

    class BruteforceOptimizer : public OptimizerInterface {
    private:
        std::size_t current_index;
    protected:
        bool UpdateSpaceWithLatestResults() override;

        SuggestionChangedParameters AlgorithmMakeSuggestion(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                                 const OptimizationParameterPoints &,
                                                                                                                 void *)> &invoke, void *data) const override;

    public:
        BruteforceOptimizer(const std::string &name, const OptimizationParameters &space, const OptimizerProperties &properties, std::size_t buffer_capacity);

        virtual ~BruteforceOptimizer();
    };

}


#endif //INMOST_BRUTEFORCE_H
