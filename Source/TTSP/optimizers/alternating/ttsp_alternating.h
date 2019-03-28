//
// Created by bvdmitri on 10.02.19.
//

#ifndef INMOST_TTSP_ALTERNATING_H
#define INMOST_TTSP_ALTERNATING_H

#include <inmost_ttsp.h>

namespace TTSP {

    enum class AlternatingDirection : std::size_t {
        RIGHT = 0,
        STAY  = 1,
        LEFT  = 2,
    };

    class AlternatingParameterHandler {
    private:
        const OptimizationParameter  &parameter;
        std::size_t                  current_index;
        mutable AlternatingDirection direction;
    public:
        explicit AlternatingParameterHandler(const OptimizationParameter &parameter);

        AlternatingParameterHandler(const AlternatingParameterHandler &other);

        const OptimizationParameter &GetParameter() const;

        AlternatingDirection GetDirection() const;

        std::size_t GetCurrentIndex() const;

        std::size_t NextIndex() const;

        void NextDirection();

        void UpdateIndex(std::size_t index);
    };

    class AlternatingOptimizer : public OptimizerInterface {
    private:
        std::size_t                              current_handler_index;
        std::vector<AlternatingParameterHandler> handlers;
    protected:
        bool UpdateSpaceWithLatestResults() override;

        OptimizationAlgorithmSuggestion AlgorithmMakeSuggestion(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                                     const OptimizationParameterPoints &,
                                                                                                                     void *)> &invoke, void *data) const override;

    public:
        AlternatingOptimizer(const std::string &name, const OptimizationParameters &space, const OptimizerProperties &properties, std::size_t buffer_capacity);

        virtual ~AlternatingOptimizer();
    };

}


#endif //INMOST_TTSP_ALTERNATING_H
