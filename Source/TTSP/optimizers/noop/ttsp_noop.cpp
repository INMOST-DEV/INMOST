//
// Created by bvdmitri on 31.01.19.
//

#include "ttsp_noop.h"

namespace TTSP {

    NoopOptimizer::NoopOptimizer(const OptimizationParameters &space, const OptimizerProperties &properties, std::size_t buffer_capacity) :
            OptimizerInterface(space, properties, buffer_capacity) {}

    OptimizationParametersSuggestion NoopOptimizer::Suggest(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                                 const OptimizationParameterPoints &,
                                                                                                                 void *)> &invoke, void *data) const {
        const OptimizationParameterEntries &entries = parameters.GetParameterEntries();

        OptimizationParameterPoints output(entries.size());

        std::transform(entries.begin(), entries.end(), output.begin(), [&](const OptimizationParametersEntry &entry) {
            return std::make_pair(entry.first.GetName(), entry.first.GetDefaultValue());
        });

        return OptimizationParametersSuggestion(entries.at(0).first, parameters.GetPoints(), parameters.GetMetrics(), output);
    }

    NoopOptimizer::~NoopOptimizer() {}

}