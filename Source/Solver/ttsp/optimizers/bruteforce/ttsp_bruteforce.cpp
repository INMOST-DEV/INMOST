//
// Created by bvdmitri on 31.01.19.
//

#include "ttsp_bruteforce.h"

namespace TTSP {


    BruteforceOptimizer::BruteforceOptimizer(const OptimizationParametersSpace &space, const OptimizerProperties &properties, std::size_t buffer_capacity) :
            OptimizerInterface(space, properties, buffer_capacity) {}

    OptimizationParametersSuggestion BruteforceOptimizer::Suggest(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                                       const OptimizationParameterPoints &,
                                                                                                                       void *)> &invoke, void *data) {

        const OptimizationParameters      &parameters = space.GetParameters();
        const OptimizationParameterPoints &before     = space.GetPoints();

        OptimizationParameterPoints output(parameters.size());

        std::transform(parameters.begin(), parameters.end(), output.begin(), [&](const OptimizationParametersEntry &entry) {
            const std::vector<double> &values = entry.first.GetValues();

            double best_metrics = -1.0;
            double best_value   = 0.0;

            std::for_each(values.begin(), values.end(), [&](double value) {
                std::cout << "[TTSP] [Bruteforce] Solving with " << entry.first.GetName() << " = " << value << "\t\t";

                const OptimizationParameterPoints &after = space.GetPointsWithChangedParameter(entry.first, value);

                OptimizationFunctionInvokeResult result = invoke(before, after, data);

                bool   is_solved = result.first;
                double metrics   = result.second;

                if (is_solved && (best_metrics < 0 || metrics < best_metrics)) {
                    best_metrics = metrics;
                    best_value   = value;
                }

                std::cout << "| Metrics = " << metrics << "\t" << is_solved << std::endl;
            });

            return std::make_pair(entry.first.GetName(), best_value);
        });

        return OptimizationParametersSuggestion(parameters.at(0).first, GetCurrentPoints(), output);
    }

    const OptimizationParameterPoints BruteforceOptimizer::GetCurrentPoints() const noexcept {
        return space.GetPoints();
    }

    BruteforceOptimizer::~BruteforceOptimizer() {}
}