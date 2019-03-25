
#include <iostream>
#include <dynamic_function.h>
#include <static_sin.h>
#include <dynamic_x2.h>

#include <inmost_ttsp.h>

using namespace TTSP;

int main(int argc, char **argv) {

    int rank = 0, size = 1;

#if defined(USE_MPI)
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Get the rank of the current process
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get the total number of processors used
#endif

    bool use_second_parameter = false;
    int  max_iterations       = 500;

    double x_default = 0.0;
    double y_default = 0.0;

    double x_step = 0.05;
    double y_step = 0.05;

    std::string optimizer_type = "noop";
    std::string function_type  = "sin1d";
    std::string out_file = "";

    //Parse argv parameters
    if (argc == 1) goto help_message;
    int i;
    for (i = 1; i < argc; i++) {
        //Print help message and exit
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            help_message:
            if (rank == 0) {
                std::cout << "Help message: " << std::endl;
                std::cout << "Command line options: " << std::endl;
                std::cout << "Required: " << std::endl;
                std::cout << "-o, --optimizer <Optimizer type>" << std::endl;
                std::cout << "Optional: " << std::endl;
                std::cout << "-f, --function <function type>" << std::endl;
                std::cout << "-usp, -y, --use-second-parameter" << std::endl;
                std::cout << "-xd, --x-default <x default value>" << std::endl;
                std::cout << "-yd, --y-default <y default value>" << std::endl;
                std::cout << "-xs, --x-step <x step>" << std::endl;
                std::cout << "-ys, --y-step <y step>" << std::endl;
                std::cout << "-mi, --max-iterations <max iterations>" << std::endl;
                std::cout << "-of, --out-file <out file path>" << std::endl;
                std::cout << "-w, --wait" << std::endl;
                std::cout << "  Available optimizers:" << std::endl;
                std::vector<std::string> availableOptimizers = TTSP::Optimizers::GetAvailableOptimizers();

                for (auto it = availableOptimizers.begin(); it != availableOptimizers.end(); ++it) {
                    std::cout << "      " << *it << std::endl;
                }
            }
            return 0;
        }
        if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--optimizer") == 0) {
            optimizer_type = std::string(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--function") == 0) {
            function_type = std::string(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-mi") == 0 || strcmp(argv[i], "--max-iterations") == 0) {
            max_iterations = std::atoi(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-y") == 0 || strcmp(argv[i], "-usp") == 0 || strcmp(argv[i], "--use-second-parameter") == 0) {
            use_second_parameter = true;
            continue;
        }
        if (strcmp(argv[i], "-xd") == 0 || strcmp(argv[i], "--x-default") == 0) {
            x_default = std::atoi(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-yd") == 0 || strcmp(argv[i], "--y-default") == 0) {
            y_default = std::atoi(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-xs") == 0 || strcmp(argv[i], "--x-step") == 0) {
            x_step = std::atoi(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-ys") == 0 || strcmp(argv[i], "--y-step") == 0) {
            y_step = std::atoi(argv[i + 1]);
            i++;
            continue;
        }
        if (strcmp(argv[i], "-of") == 0 || strcmp(argv[i], "--out-file") == 0) {
            out_file = std::string(argv[i + 1]);
            i++;
            continue;
        }
    }

    DynamicFunction *f = nullptr;
    if (function_type == "sin1d") {
        f = new StaticSin();
    } else if (function_type == "x21d") {
        f = new DynamicX2();
    } else {
        std::cout << "Available function types:" << std::endl;
        std::cout << "    sin1d - Static sin with x parameter only" << std::endl;
        std::cout << "    x21d  - Dynamic x2 with x parameter only" << std::endl;
        return 0;
    }


    OptimizationParameter x("x", f->GetXRange(), x_step, x_default);
    OptimizationParameter y("y", f->GetYRange(), y_step, y_default);

    TTSP::OptimizationParameterEntries entries;
    entries.emplace_back(std::make_pair(x, x.GetDefaultValue()));
    if (use_second_parameter) {
        entries.emplace_back(std::make_pair(y, y.GetDefaultValue()));
    }

    TTSP::OptimizationParameters parameters(entries, -1.0);

    if (!TTSP::Optimizers::IsOptimizerAvailable(optimizer_type)) {
        if (rank == 0) {
            std::cout << "Optimizer " << optimizer_type << " not found" << std::endl;
            std::cout << "  Available optimizers:" << std::endl;
            std::vector<std::string> availableOptimizers = TTSP::Optimizers::GetAvailableOptimizers();

            for (auto it = availableOptimizers.begin(); it != availableOptimizers.end(); ++it) {
                std::cout << "      " << *it << std::endl;
            }
        }
        std::exit(0);
    }

    OptimizerInterface *opt = Optimizers::GetOptimizer(optimizer_type, parameters, OptimizerProperties(), 30);

    opt->SetRestartStrategy(OptimizerRestartStrategy::RESTART_STRATEGY_WITH_BEST, 10);

    int    iteration    = 0;
    double best_metrics = -1.0;

    FILE *output = nullptr;
    if (out_file != "") {
        output = std::fopen(out_file.c_str(), "w");
    }

    while (iteration < max_iterations) {

        OptimizationParametersSuggestion suggestion = opt->Suggest();
        OptimizationParameterPoints      points     = suggestion.GetPointsAfter();

        double xv = points.at(0).GetValue();
        double yv = use_second_parameter ? points.at(1).GetValue() : 0.0;

        double metrics = f->invoke(xv, yv, iteration);

        opt->SaveResult(suggestion, metrics, true);

        if (rank == 0) {

            std::pair<double, double> min_point = f->GetMinimumPoint(iteration);

            double mx = min_point.first;
            double my = min_point.second;

            std::printf("%4d: f( %.4lf , %.4lf ) = %e | f_min(%.4lf, %.4lf) = %e\n", iteration, xv, yv, metrics, mx, my, f->GetMinimumValue(iteration));

            if (output != nullptr) {
                std::fprintf(output, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", iteration, xv, yv, metrics, mx, my, f->GetMinimumValue(iteration));
            }
        }

        iteration += 1;
    }

    if (output != nullptr) {
        std::fclose(output);
    }

    return 0;
}