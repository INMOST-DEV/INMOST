//
// Created by bvdmitri on 31.01.19.
//

#include "noop.h"

namespace INMOST {

    NoopOptimizer::NoopOptimizer(const std::string &name, const OptimizationParameters &space, const OptimizerProperties &properties, std::size_t buffer_capacity) :
            OptimizerInterface(name, space, properties, buffer_capacity) {}

    SuggestionChangedParameters NoopOptimizer::AlgorithmMakeSuggestion(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                                            const OptimizationParameterPoints &,
                                                                                                                            void *)> &invoke, void *data) const {
        return std::vector<SuggestionChangedParameter>{SuggestionChangedParameter(0, parameters.GetParameter(0).GetDefaultValue())};
    }

    NoopOptimizer::~NoopOptimizer() {}

}