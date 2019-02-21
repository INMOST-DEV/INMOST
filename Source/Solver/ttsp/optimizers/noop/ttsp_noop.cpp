//
// Created by bvdmitri on 31.01.19.
//

#include "ttsp_noop.h"

namespace TTSP {

    NoopOptimizer::NoopOptimizer(const OptimizationParametersSpace &space, const OptimizerProperties &properties, std::size_t buffer_capacity) :
            OptimizerInterface(space, properties, buffer_capacity) {}

    OptimizationParameterPoints NoopOptimizer::MakeOptimizationIteration(INMOST::Solver &solver, INMOST::Sparse::Matrix &matrix,
                                                                         INMOST::Sparse::Vector &RHS) {
        const OptimizationParameters &parameters = space.GetParameters();

        OptimizationParameterPoints output(parameters.size());

        std::transform(parameters.begin(), parameters.end(), output.begin(), [&](const OptimizationParametersEntry &entry) {
            return std::make_pair(entry.first.GetName(), entry.first.GetDefaultValue());
        });

        return output;
    }

    NoopOptimizer::~NoopOptimizer() {}

}