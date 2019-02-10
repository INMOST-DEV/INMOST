//
// Created by bvdmitri on 31.01.19.
//

#include "ttsp_bruteforce.h"

namespace TTSP {


    BruteforceOptimizer::BruteforceOptimizer(const OptimizationParametersSpace &space) : OptimizerInterface(space, 10) {}

    OptimizationParameterPoints BruteforceOptimizer::MakeOptimizationIteration(INMOST::Solver &solver, INMOST::Sparse::Matrix &matrix,
                                                                               INMOST::Sparse::Vector &RHS) {

        const OptimizationParameters &parameters = space.GetParameters();

        OptimizationParameterPoints output(parameters.size());

        std::transform(parameters.begin(), parameters.end(), output.begin(), [&](const OptimizationParametersEntry &entry) {
            const std::vector<double> &values = entry.first.GetValues();

            double tmp_time = 0.0;
            double best_time = -1.0;
            double best_value = 0.0;

            std::for_each(values.begin(), values.end(), [&](double value) {
                std::cout << "[TTSP] [Bruteforce] Solving with " << entry.first.GetName() << " = " << value << "\t\t";

                solver.SetParameter(entry.first.GetName(), INMOST::to_string(value));

                INMOST::Sparse::Vector SOL("SOL", RHS.GetFirstIndex(), RHS.GetLastIndex());
                std::fill(SOL.Begin(), SOL.End(), 0.0);

                INMOST::MPIBarrier();

                tmp_time = Timer();
                solver.SetMatrix(matrix);
                bool is_solved = solver.Solve(RHS, SOL);

                INMOST::MPIBarrier();

                double time = Timer() - tmp_time;

                if (is_solved && (best_time < 0 || time < best_time)) {
                    best_time = time;
                    best_value = value;
                }

                std::cout << "| Time = " << time << "\t" << is_solved << std::endl;
            });

            solver.SetParameter(entry.first.GetName(), INMOST::to_string(best_value));

            return std::make_pair(entry.first.GetName(), best_value);
        });

        return output;
    }

    BruteforceOptimizer::~BruteforceOptimizer() {}
}