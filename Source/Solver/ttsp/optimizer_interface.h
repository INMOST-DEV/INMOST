//
// Created by bvdmitri on 31.01.19.
//

#ifndef INMOST_OPTIMIZER_INTERFACE_H
#define INMOST_OPTIMIZER_INTERFACE_H

#include <inmost_solver.h>
#include <cstdint>
#include "optimization_parameters.h"
#include "../../Misc/utils.h"

namespace TTSP {

    class OptimizerInterface {
    protected:
        OptimizationParametersSpace space;

        void TimeBarrier() const {
#if defined(USE_MPI)
            MPI_Barrier(MPI_COMM_WORLD);
#endif
        };
    public:
        OptimizerInterface(const OptimizationParametersSpace &space) : space(space) {};

        void Solve(INMOST::Solver &solver, INMOST::Sparse::Matrix &matrix,
                   INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL) {
            const OptimizationParameterPoints &current = this->GetParametersCurrentValues();

            std::for_each(current.begin(), current.end(), [&solver](const OptimizationParameterPoint &p) {
                solver.SetParameter(p.first, INMOST::to_string(p.second));
            });

            const OptimizationParameterPoints &update = this->RequestNewParameters(current, solver, matrix, RHS, SOL);

            space.Update(update);
        };

        virtual OptimizationParameterPoints RequestNewParameters(
                const OptimizationParameterPoints &current,
                INMOST::Solver &solver, INMOST::Sparse::Matrix &matrix,
                INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL) const = 0;

        const OptimizationParameterPoints GetParametersCurrentValues() const {
            return space.GetPoints();
        }

        virtual ~OptimizerInterface() {};
    };

}

#endif // INMOST_OPTIMIZER_INTERFACE_H