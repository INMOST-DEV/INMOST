//
// Created by bvdmitri on 31.01.19.
//

#ifndef INMOST_TTSP_BRUTEFORCE_H
#define INMOST_TTSP_BRUTEFORCE_H

#include <Source/Solver/ttsp/optimizer_interface.h>

namespace TTSP {

    class BruteforceOptimizer : public OptimizerInterface {
    public:
        BruteforceOptimizer(const OptimizationParametersSpace &parameters);

        OptimizationParameterPoints RequestNewParameters(const OptimizationParameterPoints &current,
                                                         INMOST::Solver &solver, INMOST::Sparse::Matrix &matrix,
                                                         INMOST::Sparse::Vector &RHS,
                                                         INMOST::Sparse::Vector &SOL) const override;

        virtual ~BruteforceOptimizer();
    };

}


#endif //INMOST_TTSP_BRUTEFORCE_H
