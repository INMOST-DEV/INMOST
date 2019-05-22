//
// Created by bvdmitri on 31.01.19.
//

#include "bruteforce.h"

namespace INMOST {


    BruteforceOptimizer::BruteforceOptimizer(const std::string &name, const OptimizationParameters &space, const OptimizerProperties &properties, std::size_t buffer_capacity) :
            OptimizerInterface(name, space, properties, buffer_capacity), current_index(0) {}

    SuggestionChangedParameters BruteforceOptimizer::AlgorithmMakeSuggestion(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                                                  const OptimizationParameterPoints &,
                                                                                                                                  void *)> &invoke, void *data) const {
        const OptimizationParameterPoints &before    = parameters.GetPoints();
        const OptimizationParameter       &parameter = parameters.GetParameter(current_index);

        double best_metrics = -1.0;
        double best_value   = 0.0;

        std::for_each(parameter.GetValues().cbegin(), parameter.GetValues().cend(), [&](double value) {
            std::cout << "[TTSP] [Bruteforce] Solving with " << parameter.GetName() << " = " << value << "\t\t";

            const OptimizationParameterPoints &after = parameters.GetPointsWithChangedParameter(parameter, value);

            OptimizationFunctionInvokeResult result = invoke(before, after, data);

            bool   is_solved = result.first;
            double metrics   = result.second;

            if (is_solved && (best_metrics < 0 || metrics < best_metrics)) {
                best_metrics = metrics;
                best_value   = value;
            }

            std::cout << "| Metrics = " << metrics << "\t" << is_solved << std::endl;
        });

        return std::vector<SuggestionChangedParameter>{SuggestionChangedParameter(current_index, best_value)};
    }

    bool BruteforceOptimizer::UpdateSpaceWithLatestResults() {
        const OptimizationParameterResult &last = results.at(0);
        if (last.IsGood()) {
            parameters.Update(current_index, last.GetPointsAfter().at(current_index).GetValue(), last.GetMetricsAfter());
        }
        return true;
    }

    BruteforceOptimizer::~BruteforceOptimizer() {}
}