//
// Created by bvdmitri on 31.01.19.
//

#include "ttsp_bruteforce.h"

namespace TTSP {

    BruteforceOptimizer::BruteforceOptimizer(const OptimizationParametersSpace &parameters) :
            OptimizerInterface(parameters) {}


    BruteforceOptimizer::~BruteforceOptimizer() {}

    OptimizationParameterPoints BruteforceOptimizer::RequestNewParameters(const OptimizationParameterPoints &current,
                                                                          INMOST::Solver &solver,
                                                                          INMOST::Sparse::Matrix &matrix,
                                                                          INMOST::Sparse::Vector &RHS,
                                                                          INMOST::Sparse::Vector &SOL) const {

        const OptimizationParameters &parameters = space.GetParameters();

        OptimizationParameterPoints best(parameters.size());

        std::transform(parameters.begin(), parameters.end(), best.begin(), [&](const OptimizationParametersEntry &p) {
            const std::vector<double> &values = p.first.GetValues();
            double timer = 0.0;
            double best_time = -1.0;

            double best_value = 0.0;

            std::for_each(values.begin(), values.end(), [&](double value) {
                std::cout << "[TTSP] [Bruteforce] Solving with " << p.first.GetName() << " = " << value << "\t\t";

                solver.SetParameter(p.first.GetName(), INMOST::to_string(value));

                INMOST::Sparse::Vector RHSCopy(RHS);
                INMOST::Sparse::Vector SOLCopy(SOL);

                TimeBarrier();
                timer = Timer();
                solver.SetMatrix(matrix);
                bool isSolved = solver.Solve(RHSCopy, SOLCopy);
                TimeBarrier();

                double time = Timer() - timer;

                if (isSolved && (best_time < 0 || time < best_time)) {
                    best_time = time;
                    best_value = value;
                }

                std::cout << "| Time = " << time << "\t" << isSolved << std::endl;

                solver.SetParameter(p.first.GetName(), INMOST::to_string(best_value));
            });

            return std::make_pair(p.first.GetName(), best_value);
        });

        std::for_each(current.begin(), current.end(), [&solver](const OptimizationParameterPoint &point) {
            solver.SetParameter(point.first, INMOST::to_string(point.second));
        });

        solver.SetMatrix(matrix);
        solver.Solve(RHS, SOL);

        return best;
    }
}