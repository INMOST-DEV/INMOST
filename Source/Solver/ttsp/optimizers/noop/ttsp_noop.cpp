//
// Created by bvdmitri on 31.01.19.
//

#include "ttsp_noop.h"

namespace TTSP {

    NoopOptimizer::NoopOptimizer(const OptimizationParametersSpace &space, const OptimizerProperties &properties, std::size_t buffer_capacity) :
            OptimizerInterface(space, properties, buffer_capacity) {}

    OptimizationParametersSuggestion NoopOptimizer::Suggest(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                                 const OptimizationParameterPoints &,
                                                                                                                 void *)> &invoke, void *data) {
        const OptimizationParameters &parameters = space.GetParameters();

        OptimizationParameterPoints output(parameters.size());

        std::transform(parameters.begin(), parameters.end(), output.begin(), [&](const OptimizationParametersEntry &entry) {
            return std::make_pair(entry.first.GetName(), entry.first.GetDefaultValue());
        });

        return OptimizationParametersSuggestion(parameters.at(0).first, GetCurrentPoints(), output);
    }

    const OptimizationParameterPoints NoopOptimizer::GetCurrentPoints() const noexcept {
        return space.GetPoints();
    }

    NoopOptimizer::~NoopOptimizer() {}

}