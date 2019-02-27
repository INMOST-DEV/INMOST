//
// Created by bvdmitri on 31.01.19.
//

#ifndef INMOST_TTSP_DEFAULT_H
#define INMOST_TTSP_DEFAULT_H

#include <Source/Solver/ttsp/ttsp.h>

namespace TTSP {

    class NoopOptimizer : public OptimizerInterface {
    public:
        NoopOptimizer(const OptimizationParametersSpace &space, const OptimizerProperties &properties, std::size_t buffer_capacity);

        OptimizationParametersSuggestion Suggest(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                      const OptimizationParameterPoints &,
                                                                                                      void *)> &invoke, void *data) override;

        virtual ~NoopOptimizer();
    };

}

#endif //INMOST_TTSP_DEFAULT_H