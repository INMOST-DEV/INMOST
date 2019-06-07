#include "inmost.h"
#include "inmost_solver.h"
#include "inmost_optimizer.h"
#include "../../Misc/utils.h"

namespace INMOST {

    namespace TTSP {

#if defined(USE_OPTIMIZER)
        static bool                                                ttsp_disabled = true;
        static std::map<std::string, INMOST::OptimizerInterface *> optimizers    = std::map<std::string, INMOST::OptimizerInterface *>();
#endif

        void Enable() {
#if defined(USE_OPTIMIZER)
            ttsp_disabled = false;
#endif
        }

        void Disable() {
#if defined(USE_OPTIMIZER)
            ttsp_disabled = true;
#endif
        }

        bool isEnabled() {
#if defined(USE_OPTIMIZER)
            return !ttsp_disabled;
#else
            return false;
#endif
        }

        bool isDisabled() {
#if defined(USE_OPTIMIZER)
            return ttsp_disabled;
#else
            return true;
#endif
        }

#if defined(USE_OPTIMIZER)

        OptimizerInterface *GetOrCreateOptimizer(const std::string &name, const std::string &type) {
            OptimizerInterface *optimizer = nullptr;

            const std::string &key = name + "_" + type;
            if (optimizers.find(key) == optimizers.end()) {

                // TODO Here parameters for fcbiilu2 only
                OptimizationParameter        tau("tau", std::make_pair(-3, -1.0), 0.05, -2.0, OptimizationParameterType::PARAMETER_TYPE_EXPONENT);
                OptimizationParameterEntries entries;

                entries.emplace_back(std::make_pair(tau, tau.GetDefaultValue()));

                OptimizerProperties properties;
                properties["tau:use_closest"]  = "false";
                properties["tau:strict_bound"] = "false";

                OptimizationParameters parameters(entries, -1.0);

                optimizer = Optimizers::GetOptimizer(name, type, parameters, properties, 15);
                optimizers[key] = optimizer;
            } else {
                optimizer = optimizers[key];
            }
            return optimizer;
        }

#endif

        void SolverOptimize(const std::string &name, const std::string &type, Solver &solver) {
#if defined(USE_OPTIMIZER)
            if (isEnabled()) {
                OptimizerInterface *optimizer = GetOrCreateOptimizer(name, type);

                const OptimizationParametersSuggestion &suggestion = optimizer->Suggest();

                std::for_each(suggestion.GetPointsAfter().begin(), suggestion.GetPointsAfter().end(), [&solver](const OptimizationParameterPoint &point) {
                    solver.SetParameter(point.GetName(), INMOST::to_string(point.GetValue()));
                });
            }
#endif
        }

        void SolverOptimizeSaveResult(const std::string &name, const std::string &type, double metrics, bool is_good) {
#if defined(USE_OPTIMIZER)
            if (isEnabled()) {
                OptimizerInterface *optimizer = GetOrCreateOptimizer(name, type);

                auto last_suggestion = optimizer->GetLastSuggestion();
                optimizer->SaveResult(last_suggestion, metrics, is_good);
            }
#endif
        }

        void DestroySavedOptimizer(const std::string &name, const std::string &type) {
#if defined(USE_OPTIMIZER)
            if (optimizers.find(name) != optimizers.end()) {
                delete optimizers[name];
                optimizers.erase(name);
            }
#endif
        }

    }

}