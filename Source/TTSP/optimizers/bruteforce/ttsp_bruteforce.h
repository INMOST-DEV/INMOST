//
// Created by bvdmitri on 31.01.19.
//

#ifndef INMOST_TTSP_BRUTEFORCE_H
#define INMOST_TTSP_BRUTEFORCE_H

#include <inmost_ttsp.h>

namespace TTSP {

    class BruteforceOptimizer : public OptimizerInterface {
    private:
        std::size_t current_index;
    protected:
        bool UpdateSpaceWithLatestResults() override;

        OptimizationAlgorithmSuggestion AlgorithmMakeSuggestion(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                                     const OptimizationParameterPoints &,
                                                                                                                     void *)> &invoke, void *data) const override;

    public:
        BruteforceOptimizer(const std::string &name, const OptimizationParameters &space, const OptimizerProperties &properties, std::size_t buffer_capacity);

        virtual ~BruteforceOptimizer();
    };

}


#endif //INMOST_TTSP_BRUTEFORCE_H
