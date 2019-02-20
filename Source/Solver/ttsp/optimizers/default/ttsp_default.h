//
// Created by bvdmitri on 31.01.19.
//

#ifndef INMOST_TTSP_DEFAULT_H
#define INMOST_TTSP_DEFAULT_H

#include <Source/Solver/ttsp/ttsp.h>

namespace TTSP {

    class NoopOptimizer : public OptimizerInterface {
    public:
        NoopOptimizer(const OptimizationParametersSpace &space);

        OptimizationParameterPoints MakeOptimizationIteration(INMOST::Solver &solver, INMOST::Sparse::Matrix &matrix,
                                                              INMOST::Sparse::Vector &RHS) override;

        virtual ~NoopOptimizer();
    };

}

#endif //INMOST_TTSP_DEFAULT_H