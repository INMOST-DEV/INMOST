//
// Created by bvdmitri on 31.01.19.
//

#ifndef INMOST_TTSP_DEFAULT_H
#define INMOST_TTSP_DEFAULT_H

#include <inmost_optimizer.h>

namespace INMOST {

    class NoopOptimizer : public OptimizerInterface {
    protected:
        SuggestionChangedParameters AlgorithmMakeSuggestion() const override;

    public:
        NoopOptimizer(const std::string &name, const OptimizationParameters &space, const OptimizerProperties &properties, std::size_t buffer_capacity);

        virtual ~NoopOptimizer();
    };

}

#endif //INMOST_TTSP_DEFAULT_H