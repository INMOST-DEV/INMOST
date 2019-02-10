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
        const OptimizationParameter parameter;
        std::size_t direction;
        std::size_t current_index;
    public:
        explicit AlternatingParameterHandler(const OptimizationParameter &parameter);
        AlternatingParameterHandler(const AlternatingParameterHandler &other);

        const OptimizationParameter &GetParameter() const;

        std::size_t GetDirection() const;

        std::size_t GetCurrentIndex() const;

        std::size_t NextIndex();

        void NextDirection();

        void UpdateIndex(std::size_t index);
    };

    class AlternatingOptimizer : public OptimizerInterface {
    private:
        std::size_t current_handler_index;
        std::vector<AlternatingParameterHandler> handlers;
    public:
        AlternatingOptimizer(const OptimizationParametersSpace &space);

        OptimizationParameterPoints MakeOptimizationIteration(INMOST::Solver &solver, INMOST::Sparse::Matrix &matrix,
                                                              INMOST::Sparse::Vector &RHS) override;

        virtual ~AlternatingOptimizer();
    };

}



#endif //INMOST_TTSP_ALTERNATING_H
