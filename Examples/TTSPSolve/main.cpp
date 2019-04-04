#include "inmost.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>

#include <inmost_ttsp.h>
#include "series.h"
#include <Source/Misc/utils.h>

using namespace INMOST;

#if defined(USE_MPI)
#define BARRIER MPI_Barrier(MPI_COMM_WORLD)
#else
#define BARRIER
#endif


int main(int argc, char **argv) {
    int rank = 0, size = 1;

#if defined(USE_MPI)
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Get the rank of the current process
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get the total number of processors used
#endif


    {
        std::string seriesFilePath     = "";
        std::string seriesDirectory    = "";
        std::string databaseFilePath   = "";
        std::string parametersFilePath = "";
        std::string solverName         = "fcbiilu2";
        std::string optimizerType      = "bruteforce";

        bool seriesFound     = false;
        bool parametersFound = false;
        bool databaseFound   = false;
        bool waitNext        = false;

        //Parse argv parameters
        if (argc == 1) goto helpMessage;
        int i;
        for (i = 1; i < argc; i++) {
            //Print help message and exit
            if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
                helpMessage:
                if (rank == 0) {
                    std::cout << "Help message: " << std::endl;
                    std::cout << "Command line options: " << std::endl;
                    std::cout << "Required: " << std::endl;
                    std::cout << "-s, --series <Series file path>" << std::endl;
                    std::cout << "-p, --parameters <Optimizer parameters file path>" << std::endl;
                    std::cout << "Optional: " << std::endl;
                    std::cout << "-sd, --series-dir <Series directory path>" << std::endl;
                    std::cout << "-b,  --bvector <RHS vector file name>" << std::endl;
                    std::cout << "-d,  --database <Solver parameters file name>" << std::endl;
                    std::cout << "-t,  --type <Solver type name>" << std::endl;
                    std::cout << "-o,  --optt <Optimizer type name>" << std::endl;
                    std::cout << "-w,  --wait " << std::endl;
                    std::cout << "  Available solvers:" << std::endl;
                    Solver::Initialize(NULL, NULL, NULL);
                    std::vector<std::string> availableSolvers = Solver::getAvailableSolvers();
                    for (auto                it               = availableSolvers.begin(); it != availableSolvers.end(); ++it) {
                        std::cout << "      " << *it << std::endl;
                    }
                    std::cout << "  Available optimizers:" << std::endl;
                    std::vector<std::string> availableOptimizers = TTSP::Optimizers::GetAvailableOptimizers();
                    for (auto                it                  = availableOptimizers.begin(); it != availableOptimizers.end(); ++it) {
                        std::cout << "      " << *it << std::endl;
                    }
                    Solver::Finalize();
                }
                return 0;
            }
            //Series file name found with -s or --series options
            if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--series") == 0) {
                seriesFound      = true;
                seriesFilePath   = std::string(argv[i + 1]);
                FILE *seriesFile = fopen(seriesFilePath.c_str(), "r");
                if (seriesFile == NULL) {
                    if (rank == 0) {
                        std::cout << "Series file not found: " << argv[i + 1] << std::endl;
                        exit(1);
                    }
                } else {
                    if (rank == 0) {
                        std::cout << "Series file found: " << argv[i + 1] << std::endl;
                    }
                }
                fclose(seriesFile);
                i++;
                continue;
            }
            if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--parameters") == 0) {
                parametersFound      = true;
                parametersFilePath   = std::string(argv[i + 1]);
                FILE *parametersFile = fopen(parametersFilePath.c_str(), "r");
                if (parametersFile == NULL) {
                    if (rank == 0) {
                        std::cout << "Parameters file not found: " << argv[i + 1] << std::endl;
                        exit(1);
                    }
                } else {
                    if (rank == 0) {
                        std::cout << "Series file found: " << argv[i + 1] << std::endl;
                    }
                }
                fclose(parametersFile);
                i++;
                continue;
            }
            //Series directory path found with -sd or --series-dir options
            if (strcmp(argv[i], "-sd") == 0 || strcmp(argv[i], "--series-dir") == 0) {
                seriesDirectory = std::string(argv[i + 1]);
                if (rank == 0) {
                    std::cout << "Series directory prefix found: " << argv[i + 1] << std::endl;
                }
                i++;
                continue;
            }
            //Parameters file name found with -d or --database options
            if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--database") == 0) {
                if (rank == 0) {
                    std::cout << "Solver parameters file found: " << argv[i + 1] << std::endl;
                }
                databaseFound    = true;
                databaseFilePath = std::string(argv[i + 1]);
                i++;
                continue;
            }
            //Solver type found with -t ot --type options
            if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--type") == 0) {
                if (rank == 0) {
                    std::cout << "Solver type index found: " << argv[i + 1] << std::endl;
                }
                solverName = std::string(argv[i + 1]);
                i++;
                continue;
            }
            //Optimizer type found with -o ot --optt options
            if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--optt") == 0) {
                if (rank == 0) {
                    std::cout << "Optimizer type index found: " << argv[i + 1] << std::endl;
                }
                optimizerType = std::string(argv[i + 1]);
                i++;
                continue;
            }
            //Wait for each iteration
            if (strcmp(argv[i], "-w") == 0 || strcmp(argv[i], "--wait") == 0) {
                waitNext = true;
                continue;
            }
        }

        if (!seriesFound) {
            if (rank == 0) {
                std::cout <<
                          "Series file not found, you can specify series file name using -s or --series options, otherwise specify -h option to see all options, exiting...";
            }
            return -1;
        }

        if (!parametersFound) {
            if (rank == 0) {
                std::cout <<
                          "Parameters file not found, you can specify series file name using -p or --parameters options, otherwise specify -h option to see all options, exiting...";
            }
            return -1;
        }

        // Initialize series
        MatrixSeries series(seriesFilePath, seriesDirectory);

        if (series.end()) {
            if (rank == 0) {
                std::cout <<
                          "Series file found, but it looks empty or invalid, please check series file, exiting...";
            }
            return -1;
        }

        // Initialize the linear solver in accordance with args
        Solver::Initialize(&argc, &argv, databaseFound ? databaseFilePath.c_str() : NULL);

        if (!Solver::isSolverAvailable(solverName)) {
            if (rank == 0) std::cout << "Solver " << solverName << " is not available" << std::endl;
            Solver::Finalize();
            exit(1);
        }

        Solver solver = Solver(solverName, "test");

        solver.SetVerbosityLevel(SolverVerbosityLevel::Level0);
        solver.SetParameter("eps", "1e-12");

        if (rank == 0) std::cout << "Solving with " << solverName << std::endl;

        if (!TTSP::Optimizers::IsOptimizerAvailable(optimizerType)) {
            if (rank == 0) {
                std::cout << "Optimizer " << optimizerType << " not found" << std::endl;
                std::cout << "  Available optimizers:" << std::endl;
                std::vector<std::string> availableOptimizers = TTSP::Optimizers::GetAvailableOptimizers();
                for (auto                it                  = availableOptimizers.begin(); it != availableOptimizers.end(); ++it) {
                    std::cout << "      " << *it << std::endl;
                }
            }
            std::exit(0);
        }

        TTSP::OptimizersConfiguration::FromFile(parametersFilePath);

        TTSP::OptimizerProperties properties;

        properties["tau:use_closest"]  = "false";
        properties["tau:strict_bound"] = "false";

        properties["eps:use_closest"]  = "false";
        properties["eps:strict_bound"] = "false";

        TTSP::Optimizers::SaveOptimizerOrReplace("test", optimizerType, properties, 15);

        TTSP::OptimizerInterface *topt = TTSP::Optimizers::GetSavedOptimizer("test");

        topt->SetVerbosityLevel(TTSP::OptimizerVerbosityLevel::Level1);
        topt->SetRestartStrategy(TTSP::OptimizerRestartStrategy::RESTART_STRATEGY_WITH_BEST, 10);

        double metrics_total = 0.0;

        while (!series.end()) {

            std::pair<const char *, const char *> next = series.next();

            INMOST::Sparse::Matrix matrix("A");
            INMOST::Sparse::Vector rhs("b");
            INMOST::Sparse::Vector x("x");

            INMOST_DATA_ENUM_TYPE mbeg, mend;
            matrix.GetInterval(mbeg, mend);

            x.SetInterval(mbeg, mend);
            for (int k = mbeg;
                 k < mend; ++k) {
                x[k] = 0.0;
            }

            if (next.second != nullptr) {
                rhs.Load(next.second);
            } else {
                rhs.SetInterval(mbeg, mend);
                for (int k = mbeg; k < mend; ++k) rhs[k] = 1.0;
            }

            matrix.Load(next.first);

            auto invoke = [&solver, &matrix, &rhs, &x](const TTSP::OptimizationParameterPoints &before, const TTSP::OptimizationParameterPoints &after,
                                                       void *data) -> TTSP::OptimizationFunctionInvokeResult {

                std::for_each(after.begin(), after.end(), [&solver](const TTSP::OptimizationParameterPoint &point) {
                    solver.SetParameter(point.GetName(), INMOST::to_string(point.GetValue()));
                });

                INMOST::Sparse::Vector SOL("SOL", rhs.GetFirstIndex(), rhs.GetLastIndex());
                std::fill(SOL.Begin(), SOL.End(), 0.0);

                INMOST::MPIBarrier();

                double tmp_time = Timer();
                solver.SetMatrix(matrix);
                bool is_solved = solver.Solve(rhs, SOL);
                INMOST::MPIBarrier();

                double time = Timer() - tmp_time;

                std::for_each(before.begin(), before.end(), [&solver](const TTSP::OptimizationParameterPoint &point) {
                    solver.SetParameter(point.GetName(), INMOST::to_string(point.GetValue()));
                });

                return std::make_pair(is_solved, time);
            };

            TTSP::OptimizerInterface *optimizer = TTSP::Optimizers::GetSavedOptimizer("test");

            const TTSP::OptimizationParametersSuggestion &suggestion = optimizer->Suggest(invoke, nullptr);

            const TTSP::OptimizationFunctionInvokeResult &result = invoke(suggestion.GetPointsBefore(), suggestion.GetPointsAfter(), nullptr);

            bool   is_good = result.first;
            double metrics = result.second;

            metrics_total += metrics;

            optimizer->SaveResult(suggestion, metrics, is_good);

            if (rank == 0 && waitNext) {
                std::cin.get();
            }

            INMOST::MPIBarrier();
        }

        std::cout << "Metrics total from " << series.size() << " iterations: " << metrics_total << " (mean = " << metrics_total / series.size() << ")" << std::endl;
    }

    Solver::Finalize(); // Finalize solver and close MPI activity
    return 0;
}
