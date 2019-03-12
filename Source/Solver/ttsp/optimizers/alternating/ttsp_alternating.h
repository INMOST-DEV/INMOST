//
// Created by bvdmitri on 10.02.19.
//

#ifndef INMOST_TTSP_ALTERNATING_H
#define INMOST_TTSP_ALTERNATING_H

#include <Source/Solver/ttsp/ttsp.h>

namespace TTSP {

    enum AlternatingDirection : std::size_t {
        RIGHT = 0,
        STAY1 = 1,
        LEFT  = 2,
        STAY2 = 3,
    };

    class AlternatingParameterHandler {
    private:
        const OptimizationParameter &parameter;
        std::size_t                 current_index;
        mutable std::size_t         direction;
    public:
        explicit AlternatingParameterHandler(const OptimizationParameter &parameter);

        AlternatingParameterHandler(const AlternatingParameterHandler &other);

        const OptimizationParameter &GetParameter() const;

        std::size_t GetDirection() const;

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
        void UpdateSpaceWithLatestResults() override;

    public:
        AlternatingOptimizer(const OptimizationParameters &space, const OptimizerProperties &properties, std::size_t buffer_capacity);

        OptimizationParametersSuggestion Suggest(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                      const OptimizationParameterPoints &,
                                                                                                      void *)> &invoke, void *data) const override;

        virtual ~AlternatingOptimizer();
    };

}


#endif //INMOST_TTSP_ALTERNATING_H
