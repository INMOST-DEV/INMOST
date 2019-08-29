#include "inmost.h"
#include "inmost_solver.h"
#include "inmost_optimizer.h"
#include "../../Misc/utils.h"
#include "ttsp_configuration.h"

namespace INMOST {

    namespace TTSP {

#if defined(USE_OPTIMIZER)
        struct TTSPOptions {
            bool enabled;
        };

        static TTSPConfiguration                                   *configuration = nullptr;
        static std::map<std::string, TTSPOptions>                  g_options      = std::map<std::string, TTSPOptions>();
        static std::map<std::string, INMOST::OptimizerInterface *> g_optimizers   = std::map<std::string, INMOST::OptimizerInterface *>();
#endif

        void Initialize(const std::string &path) {
#if defined(USE_OPTIMIZER)
            Deinitialize();

            INMOST::TTSP::configuration = new TTSPConfiguration(path);

            std::for_each(INMOST::TTSP::configuration->GetSolvers().cbegin(), INMOST::TTSP::configuration->GetSolvers().cend(), [](const TTSPConfigurationSolverEntry &s) {
                if (s.IsEnabled()) {
                    std::for_each(s.GetPrefixes().cbegin(), s.GetPrefixes().cend(), [&s](const TTSPConfigurationSolverPrefixEntry &p) {
                        if (p.IsEnabled()) {
                            Enable(s.GetSolver(), p.GetPrefix());
                        }
                    });
                }
            });
#endif
        }

        void Deinitialize() {
#if defined(USE_OPTIMIZER)
            delete INMOST::TTSP::configuration;
            INMOST::TTSP::configuration = nullptr;
            g_options.clear();
            g_optimizers.clear();
#endif
        }

        void Enable(const std::string &solver_name, const std::string &solver_prefix) {
#if defined(USE_OPTIMIZER)
            const std::string &key = solver_name + ":" + solver_prefix;
            if (g_options.find(key) == g_options.end()) {
                TTSPOptions n_options;
                n_options.enabled = true;
                g_options[key] = n_options;
            } else {
                g_options[key].enabled = true;
            }
#endif
        }

        void Disable(const std::string &solver_name, const std::string &solver_prefix) {
#if defined(USE_OPTIMIZER)
            const std::string &key = solver_name + ":" + solver_prefix;
            if (g_options.find(key) != g_options.end()) {
                g_options[key].enabled = false;
            }
#endif
        }

        bool isEnabled(const std::string &solver_name, const std::string &solver_prefix) {
#if defined(USE_OPTIMIZER)
            const std::string &key = solver_name + ":" + solver_prefix;
            if (g_options.find(key) != g_options.end()) {
                return g_options[key].enabled;
            } else {
                return false;
            }
#else
            return false;
#endif
        }

        bool isDisabled(const std::string &solver_name, const std::string &solver_prefix) {
#if defined(USE_OPTIMIZER)
            const std::string &key = solver_name + ":" + solver_prefix;
            if (g_options.find(key) != g_options.end()) {
                return !g_options[key].enabled;
            } else {
                return true;
            }
#else
            return true;
#endif
        }

#if defined(USE_OPTIMIZER)

        OptimizerInterface *GetOrCreateOptimizer(Solver &solver) {
            OptimizerInterface *optimizer = nullptr;
            const std::string  &key       = solver.SolverName() + ":" + solver.SolverPrefix();
            if (g_optimizers.find(key) == g_optimizers.end()) {

                OptimizationParameterEntries entries;

                const TTSPConfigurationSolverPrefixEntry &entry = INMOST::TTSP::configuration->FindForSolverAndPrefix(solver.SolverName(), solver.SolverPrefix());

                std::for_each(entry.GetParameters().cbegin(), entry.GetParameters().cend(), [&entries](const TTSPConfigurationParameterEntry &e) {
                    OptimizationParameter o(e.GetName(), e.GetValues(), e.GetInitial(), e.GetParameterType());
                    entries.emplace_back(std::make_pair(o, o.GetDefaultValue()));
                });

                OptimizationParameters optparams(entries, -1.0);

                OptimizerProperties properties;

                optimizer = Optimizers::GetOptimizer(key, entry.GetOptimizer(), optparams, properties, entry.GetBufferCapacity());

                optimizer->SetVerbosityLevel(entry.GetVerbosityLevel());
                optimizer->SetRestartStrategy(OptimizerRestartStrategy::RESTART_STRATEGY_WITH_BEST, entry.GetBufferCapacity() - 1);

                g_optimizers[key] = optimizer;
            } else {
                optimizer = g_optimizers[key];
            }
            return optimizer;
        }

#endif

        void SolverOptimize(Solver &solver, bool use_last_suggestion) {
#if defined(USE_OPTIMIZER)
            if (isEnabled(solver.SolverName(), solver.SolverPrefix())) {
                OptimizerInterface *optimizer = GetOrCreateOptimizer(solver);

                if (optimizer != nullptr) {
                    const OptimizationParametersSuggestion &suggestion = use_last_suggestion ? optimizer->GetLastSuggestion() : optimizer->Suggest();

                    std::for_each(suggestion.GetPointsAfter().begin(), suggestion.GetPointsAfter().end(), [&solver](const OptimizationParameterPoint &point) {
                        solver.SetParameter(point.GetName(), INMOST::to_string(point.GetValue()));
                    });
                }
            }
#endif
        }

        void SolverOptimizeSaveResult(Solver &solver, double metrics, bool is_good) {
#if defined(USE_OPTIMIZER)
            if (isEnabled(solver.SolverName(), solver.SolverPrefix())) {
                OptimizerInterface *optimizer = GetOrCreateOptimizer(solver);
                if (optimizer != nullptr) {
                    auto last_suggestion = optimizer->GetLastSuggestion();
                    optimizer->SaveResult(last_suggestion, metrics, is_good);
                }
            }
#endif
        }

        void DestroySavedOptimizer(Solver &solver) {
#if defined(USE_OPTIMIZER)
            const std::string &key = solver.SolverName() + ":" + solver.SolverPrefix();
            if (g_optimizers.find(key) != g_optimizers.end()) {
                delete g_optimizers[key];
                g_optimizers.erase(key);
            }
#endif
        }

    }

}