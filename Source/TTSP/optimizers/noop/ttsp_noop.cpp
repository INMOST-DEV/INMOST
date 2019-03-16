//
// Created by bvdmitri on 31.01.19.
//

#include "ttsp_noop.h"

namespace TTSP {

    NoopOptimizer::NoopOptimizer(const OptimizationParameters &space, const OptimizerProperties &properties, std::size_t buffer_capacity) :
            OptimizerInterface(space, properties, buffer_capacity) {}

    OptimizationAlgorithmSuggestion NoopOptimizer::AlgorithmMakeSuggestion(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                                                const OptimizationParameterPoints &,
                                                                                                                                void *)> &invoke, void *data) const {
        return std::make_pair(0, parameters.GetParameter(0).GetDefaultValue());
    }

    NoopOptimizer::~NoopOptimizer() {}

}