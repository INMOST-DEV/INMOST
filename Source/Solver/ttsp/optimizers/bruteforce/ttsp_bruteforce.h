//
// Created by bvdmitri on 31.01.19.
//

#ifndef INMOST_TTSP_BRUTEFORCE_H
#define INMOST_TTSP_BRUTEFORCE_H

#include <Source/Solver/ttsp/ttsp.h>

namespace TTSP {

    class BruteforceOptimizer : public OptimizerInterface {
    public:
        BruteforceOptimizer(const OptimizationParametersSpace &space, const OptimizerProperties &properties, std::size_t buffer_capacity);

        OptimizationParameterPoints MakeOptimizationIteration(INMOST::Solver &solver, INMOST::Sparse::Matrix &matrix,
                                                              INMOST::Sparse::Vector &RHS) override;

        virtual ~BruteforceOptimizer();
    };

}


#endif //INMOST_TTSP_BRUTEFORCE_H
