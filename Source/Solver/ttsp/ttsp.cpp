#include "inmost.h"
#include "inmost_solver.h"
#include "inmost_optimizer.h"
#include "../../Misc/utils.h"

namespace INMOST {

    namespace TTSP {

#if defined(USE_OPTIMIZER)
        struct TTSPParameters {
            bool        disabled;
            std::string type;
            // TODO add array of specific parameters here, like tau and so on
        };

        static std::map<std::string, TTSPParameters>               g_parameters = std::map<std::string, TTSPParameters>();
        static std::map<std::string, INMOST::OptimizerInterface *> g_optimizers = std::map<std::string, INMOST::OptimizerInterface *>();
#endif

        void Enable(const std::string &name, const std::string &type) {
#if defined(USE_OPTIMIZER)
            if (g_parameters.find(name) == g_parameters.end()) {
                TTSPParameters p;

                p.disabled = false;
                p.type     = type;

                g_parameters[name] = p;
            } else {
                g_parameters[name].disabled = false;
            }
#endif
        }

        void Disable(const std::string &name, const std::string &type) {
#if defined(USE_OPTIMIZER)
            if (g_parameters.find(name) != g_parameters.end()) {
                g_parameters[name].disabled = true;
            }
#endif
        }

        bool isEnabled(const std::string &name) {
#if defined(USE_OPTIMIZER)
            if (g_parameters.find(name) != g_parameters.end()) {
                return !g_parameters[name].disabled;
            } else {
                return false;
            }
#else
            return false;
#endif
        }

        bool isDisabled(const std::string &name) {
#if defined(USE_OPTIMIZER)
            if (g_parameters.find(name) != g_parameters.end()) {
                return g_parameters[name].disabled;
            } else {
                return true;
            }
#else
            return true;
#endif
        }

#if defined(USE_OPTIMIZER)

        OptimizerInterface *GetOrCreateOptimizer(const std::string &name) {
            OptimizerInterface *optimizer = nullptr;

            if (g_optimizers.find(name) == g_optimizers.end()) {

                // TODO Here parameters for fcbiilu2 only
                OptimizationParameter        tau("tau", std::make_pair(-3, -1.0), 0.05, -2.0, OptimizationParameterType::PARAMETER_TYPE_EXPONENT);
                OptimizationParameterEntries entries;

                entries.emplace_back(std::make_pair(tau, tau.GetDefaultValue()));

                OptimizerProperties properties;
                properties["tau:use_closest"]  = "false";
                properties["tau:strict_bound"] = "false";

                OptimizationParameters optparams(entries, -1.0);

                optimizer = Optimizers::GetOptimizer(name, g_parameters[name].type, optparams, properties, 15);

                optimizer->SetVerbosityLevel(OptimizerVerbosityLevel::Level3);
                optimizer->SetRestartStrategy(OptimizerRestartStrategy::RESTART_STRATEGY_WITH_BEST, 15);

                g_optimizers[name] = optimizer;
            } else {
                optimizer = g_optimizers[name];
            }
            return optimizer;
        }

#endif

        void SolverOptimize(const std::string &name, Solver &solver, bool use_last_suggestion) {
#if defined(USE_OPTIMIZER)
            if (isEnabled(name)) {
                OptimizerInterface *optimizer = GetOrCreateOptimizer(name);

                const OptimizationParametersSuggestion &suggestion = use_last_suggestion ? optimizer->GetLastSuggestion() : optimizer->Suggest();

                std::for_each(suggestion.GetPointsAfter().begin(), suggestion.GetPointsAfter().end(), [&solver](const OptimizationParameterPoint &point) {
                    solver.SetParameter(point.GetName(), INMOST::to_string(point.GetValue()));
                });
            }
#endif
        }

        void SolverOptimizeSaveResult(const std::string &name, double metrics, bool is_good) {
#if defined(USE_OPTIMIZER)
            if (isEnabled(name)) {
                OptimizerInterface *optimizer = GetOrCreateOptimizer(name);

                auto last_suggestion = optimizer->GetLastSuggestion();
                optimizer->SaveResult(last_suggestion, metrics, is_good);
            }
#endif
        }

        void DestroySavedOptimizer(const std::string &name) {
#if defined(USE_OPTIMIZER)
            if (g_optimizers.find(name) != g_optimizers.end()) {
                delete g_optimizers[name];
                g_optimizers.erase(name);
            }
#endif
        }

    }

}