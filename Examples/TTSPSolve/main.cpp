#include "inmost.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>

#include "Source/Solver/ttsp/ttsp.h"
#include "Source/Solver/ttsp/optimizers/bruteforce/ttsp_bruteforce.h"

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
        std::string matrixFileName = "";
        std::string vectorBFileName = "";
        std::string parametersFileName = "";
        std::string solverName = "fcbiilu2";

        bool matrixFound = false;
        bool vectorBFound = false;
        bool parametersFound = false;
        bool typeFound = false;

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
                    std::cout << "-m, --matrix <Matrix file name>" << std::endl;
                    std::cout << "Optional: " << std::endl;
                    std::cout << "-b, --bvector <RHS vector file name>" << std::endl;
                    std::cout << "-d, --database <Solver parameters file name>" << std::endl;
                    std::cout << "-t, --type <Solver type name>" << std::endl;
                    std::cout << "  Available solvers:" << std::endl;
                    Solver::Initialize(NULL, NULL, NULL);
                    std::vector<std::string> availableSolvers = Solver::getAvailableSolvers();
                    for (auto it = availableSolvers.begin(); it != availableSolvers.end(); it++) {
                        std::cout << "      " << *it << std::endl;
                    }
                    Solver::Finalize();
                }
                return 0;
            }
            //Matrix file name found with -m or --matrix options
            if (strcmp(argv[i], "-m") == 0 || strcmp(argv[i], "--matrix") == 0) {
                matrixFound = true;
                matrixFileName = std::string(argv[i + 1]);
                FILE *matrixFile = fopen(matrixFileName.c_str(), "r");
                if (matrixFile == NULL) {
                    if (rank == 0) {
                        std::cout << "Matrix file not found: " << argv[i + 1] << std::endl;
                        exit(1);
                    }
                } else {
                    if (rank == 0) {
                        std::cout << "Matrix file found: " << argv[i + 1] << std::endl;
                    }
                }
                fclose(matrixFile);
                i++;
                continue;
            }
            //B vector file name found with -b or --bvector options
            if (strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--bvector") == 0) {
                if (rank == 0) {
                    std::cout << "B vector file found: " << argv[i + 1] << std::endl;
                }
                vectorBFound = true;
                vectorBFileName = std::string(argv[i + 1]);
                i++;
                continue;
            }
            //Parameters file name found with -d or --database options
            if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--database") == 0) {
                if (rank == 0) {
                    std::cout << "Solver parameters file found: " << argv[i + 1] << std::endl;
                }
                parametersFound = true;
                parametersFileName = std::string(argv[i + 1]);
                i++;
                continue;
            }
            //Solver type found with -t ot --type options
            if (strcmp(argv[i], "-t") == 0 || strcmp(argv[i], "--type") == 0) {
                if (rank == 0) {
                    std::cout << "Solver type index found: " << argv[i + 1] << std::endl;
                }
                typeFound = true;
                solverName = std::string(argv[i + 1]);
                i++;
                continue;
            }
        }

        if (!matrixFound) {
            if (rank == 0) {
                std::cout <<
                          "Matrix not found, you can specify matrix file name using -m or --matrix options, otherwise specify -h option to see all options, exiting...";
            }
            return -1;
        }

        if (!typeFound) {
            if (rank == 0) {
                std::cout <<
                          "Solver type not found in command line, you can specify solver type with -t or --type option, using INNER_ILU2 solver by default."
                          <<
                          std::endl;
            }
        }

        if (!vectorBFound) {
            if (rank == 0) {
                std::cout <<
                          "B vector not found, you can specify b vector file name with -b or --bvector option, using identity vector by default."
                          <<
                          std::endl;
            }
        }

        // Initialize the linear solver in accordance with args
        Solver::Initialize(&argc, &argv, parametersFound ? parametersFileName.c_str() : NULL);

        if (!Solver::isSolverAvailable(solverName)) {
            if (rank == 0) std::cout << "Solver " << solverName << " is not available" << std::endl;
            Solver::Finalize();
            exit(1);
        }

        Solver solver = Solver(solverName, "test");

        if (rank == 0) std::cout << "Solving with " << solverName << std::endl;

        Sparse::Matrix mat("A"); // Declare the matrix of the linear system to be solved
        Sparse::Vector b("rhs"); // Declare the right-hand side vector
        Sparse::Vector x("sol"); // Declare the solution vector

        double timer = Timer();
        mat.Load(matrixFileName); //ifinterval parameters not set, matrix will be divided automatically

        if (rank == 0) std::cout << "Load matrix time:    " << Timer() - timer << std::endl;

        timer = Timer();
        if (vectorBFound) {
            b.Load(vectorBFileName); //if interval parameters not set, matrix will be divided automatically
        } else { // Set local RHS to 1 if it was not specified
            INMOST_DATA_ENUM_TYPE mbeg, mend, k;
            mat.GetInterval(mbeg, mend);
            b.SetInterval(mbeg, mend);
            for (k = mbeg; k < mend; ++k) b[k] = 1.0;
        }

        if (rank == 0) std::cout << "Load vector time:    " << Timer() - timer << std::endl;

        TTSP::OptimizationParameter tau("tau", {1e-3, 3e-3, 5e-3, 7e-3, 1e-2, 3e-2, 5e-2, 7e-2}, 1e-3);


        TTSP::OptimizationParameters parameters;
        parameters.push_back(std::make_pair(tau, 1e-3));

        TTSP::OptimizationParametersSpace space(solverName, "test", parameters);
        TTSP::BruteforceOptimizer optimizer(space);

        //        BARRIER;
        //        timer = Timer();
        //        solver.SetMatrix(mat);       // Compute the preconditioner for the original matrix
        //        BARRIER;
        //
        //        if (rank == 0) std::cout << "Preconditioner time: " << Timer() - timer << std::endl;
        //
        //        BARRIER;
        //        timer = Timer();
        //        bool isSuccess = solver.Solve(b, x); // Solve the linear system with the previously computted preconditioner
        //        BARRIER;

        int test = 0;

        while (test < 15) {

            optimizer.Solve(solver, mat, b, x);

            std::cout << std::endl << "Best optimization parameters found for current iteration:" << std::endl;
            const TTSP::OptimizationParameterPoints &best = optimizer.GetSpace().GetPoints();
            std::for_each(best.begin(), best.end(), [](const TTSP::OptimizationParameterPoint &p) {
                std::cout << "\t" << p.first << " = " << p.second << std::endl;
            });

            std::cout << std::endl << "Optimization results buffer output:" << std::endl;
            const TTSP::OptimizationParameterResultsBuffer &results = optimizer.GetResults();

            int index = 1;
            std::for_each(results.begin(), results.end(), [&index](const TTSP::OptimizationParameterResult &result) {
                std::cout << "\t" << index++ << "\t" << " [";

                const TTSP::OptimizationParameterPoints &points = result.GetPoints();
                std::for_each(points.begin(), points.end(), [](const TTSP::OptimizationParameterPoint &point) {
                    std::cout << " " << point.first << "=" << point.second << " ";
                });

                std::cout << "] " << result.GetPreconditionerTime() << "\t" << result.GetSolveTime() << "\t" << result.GetTime() << std::endl;
            });

            if (rank == 0) {
                std::cout << std::endl
                          << "Solved with " << solver.SolverName()
                          << " on " << solver.Iterations()
                          << " iterations and " << solver.Residual()
                          << " residual. Reason: " << solver.ReturnReason()
                          << std::endl;
            }

            test += 1;
        }
    }

    Solver::Finalize(); // Finalize solver and close MPI activity
    return 0;
}
