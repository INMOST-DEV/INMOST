#include "inmost.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>

#include "Source/Solver/ttsp/ttsp.h"
#include "series.h"

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
        std::string seriesFileName     = "";
        std::string seriesDirectory    = "";
        std::string parametersFileName = "";
        std::string solverName         = "fcbiilu2";
        std::string optimizerType      = "bruteforce";

        bool seriesFound     = false;
        bool parametersFound = false;
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
                    std::cout << "-s, --series <Series file name>" << std::endl;
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
                    std::vector<std::string> availableOptimizers = TTSP::OptimizerInterface::GetAvailableOptimizers();
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
                seriesFileName   = std::string(argv[i + 1]);
                FILE *seriesFile = fopen(seriesFileName.c_str(), "r");
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
                parametersFound    = true;
                parametersFileName = std::string(argv[i + 1]);
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

        // Initialize series
        MatrixSeries series(seriesFileName, seriesDirectory);

        if (series.end()) {
            if (rank == 0) {
                std::cout <<
                          "Series file found, but it looks empty or invalid, please check series file, exiting...";
            }
            return -1;
        }

        // Initialize the linear solver in accordance with args
        Solver::Initialize(&argc, &argv, parametersFound ? parametersFileName.c_str() : NULL);

        if (!Solver::isSolverAvailable(solverName)) {
            if (rank == 0) std::cout << "Solver " << solverName << " is not available" << std::endl;
            Solver::Finalize();
            exit(1);
        }

        Solver solver = Solver(solverName, "test");

        solver.SetVerbosityLevel(SolverVerbosityLevel::Level0);

        if (rank == 0) std::cout << "Solving with " << solverName << std::endl;

        TTSP::OptimizationParameter  tau("tau", {3e-2, 5e-2, 6e-2, 7e-2, 8e-2, 9e-2, 1e-1, 2e-1, 3e-1, 5e-1, 7e-1, 9e-1}, 1e-3);
        TTSP::OptimizationParameters parameters;
        parameters.push_back(std::make_pair(tau, 1e-3));

        TTSP::OptimizerProperties properties;

        properties["tau:use_closest"]  = "false";
        properties["tau:strict_bound"] = "false";

        TTSP::OptimizationParametersSpace space(solverName, "test", parameters);

        TTSP::OptimizerInterface *optimizer = TTSP::OptimizerInterface::GetOptimizer(optimizerType, space, properties, 10);
        if (optimizer == nullptr) {
            if (rank == 0) {
                std::cout << "Optimizer " << optimizerType << " not found" << std::endl;
                std::cout << "  Available optimizers:" << std::endl;
                std::vector<std::string> availableOptimizers = TTSP::OptimizerInterface::GetAvailableOptimizers();
                for (auto                it                  = availableOptimizers.begin(); it != availableOptimizers.end(); ++it) {
                    std::cout << "      " << *it << std::endl;
                }
            }
            std::exit(0);
        }

        optimizer->SetVerbosityLevel(TTSP::OptimizerVerbosityLevel::Level1);

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
                    solver.SetParameter(point.first, INMOST::to_string(point.second));
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
                    solver.SetParameter(point.first, INMOST::to_string(point.second));
                });

                return std::make_pair(is_solved, time);
            };

            const TTSP::OptimizationParameterPoints      &before     = optimizer->GetSpace().GetPoints();
            const TTSP::OptimizationParametersSuggestion &suggestion = optimizer->Suggest(invoke, nullptr);

            const TTSP::OptimizationFunctionInvokeResult &result = invoke(before, suggestion.second, nullptr);

            bool   is_good = result.first;
            double metrics = result.second;

            optimizer->SaveResult(suggestion.first, before, suggestion.second, metrics, is_good);

            if (is_good) {
                optimizer->UpdateSpacePoints(suggestion.second);
            }

            TTSP::OptimizerVerbosityLevel verbosity = TTSP::OptimizerVerbosityLevel::Level3;

            // On Level1 print some metadata information about solution and used parameters
            if (rank == 0 && verbosity > TTSP::OptimizerVerbosityLevel::Level0) {
                std::string metadata = solver.SolutionMetadataLine("\t");
                std::for_each(suggestion.second.begin(), suggestion.second.end(), [&metadata](const TTSP::OptimizationParameterPoint &p) {
                    metadata += ("\t" + INMOST::to_string(p.second));
                });
                std::cout << metadata << std::endl;
            }

            // On Level2 also print information about next parameters
            if (rank == 0 && verbosity > TTSP::OptimizerVerbosityLevel::Level1) {
                std::cout << std::endl << "Next optimization parameters found for current iteration:" << std::endl;
                const TTSP::OptimizationParameterPoints &points = optimizer->GetSpace().GetPoints();
                std::for_each(points.begin(), points.end(), [](const TTSP::OptimizationParameterPoint &p) {
                    std::cout << "\t" << p.first << " = " << p.second << std::endl;
                });
            }

            // On Level3 also print additional information about buffer
            if (rank == 0 && verbosity > TTSP::OptimizerVerbosityLevel::Level2) {
                std::cout << std::endl << "Optimization results buffer output:" << std::endl;
                const TTSP::OptimizationParameterResultsBuffer &results = optimizer->GetResults();

                int index = 1;
                std::for_each(results.begin(), results.end(), [&index](const TTSP::OptimizationParameterResult &result) {
                    std::cout << "\t" << index++ << "\t" << " [";

                    const TTSP::OptimizationParameterPoints &before = result.GetPointsBefore();
                    std::for_each(before.begin(), before.end(), [](const TTSP::OptimizationParameterPoint &point) {
                        std::cout << " " << point.first << "=" << point.second << " ";
                    });

                    std::cout << "] >> [";

                    const TTSP::OptimizationParameterPoints &after = result.GetPointsAfter();
                    std::for_each(after.begin(), after.end(), [](const TTSP::OptimizationParameterPoint &point) {
                        std::cout << " " << point.first << "=" << point.second << " ";
                    });

                    std::cout << "] " << result.GetMetrics() << std::endl;
                });
            }

            if (rank == 0 && waitNext) {
                std::cin.get();
            }

            INMOST::MPIBarrier();
        }
    }

    Solver::Finalize(); // Finalize solver and close MPI activity
    return 0;
}
