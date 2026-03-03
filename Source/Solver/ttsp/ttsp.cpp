#include "inmost.h"
#include "inmost_solver.h"
#include "inmost_optimizer.h"
#include "../../Misc/utils.h"
#include "ttsp_configuration.h"
#if defined(USE_SOLVER)
namespace INMOST {

    namespace TTSP {

//#define USE_SOLVER_STAT
#if defined(USE_SOLVER_STAT)

        class TTSP3Entry {
        public:
            std::string solvername;
            std::string prefix;
            double      parameter;
            double      tprec;
            double      titer;
            int         count;
            int         nbad;
        
            TTSP3Entry() {solvername=""; prefix=""; parameter=0.0; tprec=0.0; titer=0.0; count=0; nbad=0;}
            TTSP3Entry(std::string prefix) {solvername=""; this->prefix=prefix; parameter=0.0; tprec=0.0; titer=0.0; count=0; nbad=0;}
            ~TTSP3Entry() {};
        };

        static TTSP3Entry e[3] = {
            INMOST::TTSP::TTSP3Entry("flow"),
            INMOST::TTSP::TTSP3Entry("tran"),
            INMOST::TTSP::TTSP3Entry("heat")
        };

        void SolverOptimizeSaveResult3(Solver &solver, double metrics, bool is_good) {
            (void)solver,(void)metrics,(void)is_good;
            int i = -1;
            std::string prefix = solver.SolverPrefix();
            if (prefix == e[0].prefix)
                i = 0;
            else if (prefix == e[1].prefix)
                i = 1;
            else if (prefix == e[2].prefix)
                i = 2;
            if (i >= 0) {
                if (e[i].count == 0) {
                    e[i].solvername = solver.SolverName();
                    e[i].parameter = 0.0; //TODO: save parameter value for further analisys
                }
                if (is_good) {
                    e[i].tprec += solver.PreconditionerTime();
                    e[i].titer += solver.IterationsTime();
                } else {
                    e[i].nbad++;
                }
                e[i].count++;
            }
        }

        static std::string to_string2(double value)
        {
            std::stringstream ss;
            ss.precision(2);
            ss << std::fixed << value;
            return ss.str();
        }

        std::string SolutionMetadataLine3() {
            bool eng = false; // English or Russian output
            bool parallel = false; // parallel or serial run
#if defined (USE_MPI)
            int nproc;
            MPI_Comm_size(MPI_COMM_WORLD, &nproc);
            if (nproc > 1) parallel = true;
#endif
            std::string str = "";
            for (int i = 0; i < 3; i++) {
                if (e[i].count) {
                    double total = e[i].tprec + e[i].titer + 1e-9;
                    double tprec = e[i].tprec / total * 100.0;
                    double titer = e[i].titer / total * 100.0;
                    str += INMOST::to_string(e[i].prefix) + "(" +
                           INMOST::to_string(e[i].solvername) + "): " +
                           to_string2(total) + (eng ? " sec. = " : " сек. = ") +
                           to_string2(tprec) + "%(prec) + " +
                           to_string2(titer) + "%(iter) " + (eng ? "for " : "для ") +
                           INMOST::to_string(e[i].count) + (eng ? " calls (" : " вызовов (") +
                           to_string2(100.0-e[i].nbad*100.0/e[i].count) + (eng ? "% successfully" : "% успешно") +
                           ")\n";
                    int stronger = 0; // 2, 1, 0, -1, -2 :: [0 : 20 : 40-50-60 : 80 : 100]
                    double t20 = 20, t40 = 40, t60 = 60, t80 = 80;
                    if (e[i].solvername == "petsc") { // preconditioner for petsc may be weaker
                        t20 = 10, t40 = 30, t60 = 50, t80 = 70;
                    } else if (e[i].solvername != "biilu2") { // preconditioner for inner INMOST solvers may be stronger
                        t20 = 25, t40 = 45, t60 = 65, t80 = 85;
                    }
                    if (0.01 <= tprec && tprec < t20)
                        stronger = 2;
                    else if (t20 <= tprec && tprec < t40)
                        stronger = 1;
                    else if (t60 < tprec && tprec <= t80)
                        stronger = -1;
                    else if (t80 < tprec)
                        stronger = -2;
                    if (e[i].nbad == 0 && stronger == 0) {
                        str += (eng) ? "[Parameters are close to optimal" :
                                       "[Параметры близки к оптимальным";
                    } else if (e[i].nbad != 0) {
                        str += (eng) ? "[Parameters can be made a little stronger due to divergence" :
                                       "[Параметры можно немного усилить из-за случаев несходимости";
                        stronger = 1;
                    } else if (stronger != 0) {
                        if (stronger == 2)
                            str += (eng) ? "[Parameters can be made stronger" :
                                           "[Параметры можно усилить";
                        else if (stronger == 1)
                            str += (eng) ? "[Parameters can be made a little stronger" :
                                           "[Параметры можно немного усилить";
                        else if (stronger == -1)
                            str += (eng) ? "[Parameters can be made a little weaker" :
                                           "[Параметры можно сделать немного слабее";
                        else if (stronger == -2)
                            str += (eng) ? "[Parameters can be made weaker" :
                                           "[Параметры можно сделать слабее";
                    }
// petsc :: sub_pc_factor_levels, pc_asm_overlap
// biilu2 :: tau kovl
// inner_ilu2 inner_ddpqiluc2 inner_mptiluc inner_mlmptiluc inner_mptilu2 :: drop_tolerance+reuse_tolerance, schwartz_overlap
                    if (e[i].solvername == "petsc" && (stronger > 0 || e[i].nbad != 0)) {
                        str += (eng ? ": increase sub_pc_factor_levels" :
                                       ": увеличить sub_pc_factor_levels"); 
                        if (parallel) str += (eng) ? " and/or pc_asm_overlap" :
                                                     " и/или pc_asm_overlap";
                    } else if (e[i].solvername == "petsc" && stronger < 0 && e[i].nbad != 0) {
                        str += (eng) ? ": reduce sub_pc_factor_levels" :
                                       ": уменьшить sub_pc_factor_levels";
                        if (parallel) str += (eng) ? " and/or pc_asm_overlap" :
                                                     " и/или pc_asm_overlap";
                    } else if (e[i].solvername == "biilu2" && stronger > 0) {
                        str += (eng) ? ": reduce tau" :
                                       ": уменьшить tau";
                        if (parallel) str += (eng) ? " and/or increase kovl" :
                                                     " и/или увеличить kovl";
                    } else if (e[i].solvername == "biilu2" && stronger < 0) {
                        str += (eng) ? ": increase tau" :
                                       ": увеличить tau";
                        if (parallel) str += (eng) ? " and/or reduce kovl" :
                                                     " и/или уменьшить kovl";
                    } else if (stronger > 0 && (e[i].solvername == "inner" || e[i].solvername == "inner_ilu2" || e[i].solvername == "inner_ddpqiluc2" || e[i].solvername == "inner_mptiluc" || e[i].solvername == "inner_mlmptiluc" || e[i].solvername == "inner_mptilu2")) {
                        str += (eng) ? ": reduce drop_tolerance and reuse_tolerance" :
                                       ": уменьшить drop_tolerance и reuse_tolerance";
                        if (parallel) str += (eng) ? " and/or increase schwartz_overlap" :
                                                     " и/или увеличить schwartz_overlap";
                    } else if (stronger < 0 && (e[i].solvername == "inner" || e[i].solvername == "inner_ilu2" || e[i].solvername == "inner_ddpqiluc2" || e[i].solvername == "inner_mptiluc" || e[i].solvername == "inner_mlmptiluc" || e[i].solvername == "inner_mptilu2")) {
                        str += (eng) ? ": increase drop_tolerance and reuse_tolerance" :
                                       ": увеличить drop_tolerance и reuse_tolerance";
                        if (parallel) str += (eng) ? " and/or reduce schwartz_overlap" :
                                                     " и/или уменьшить schwartz_overlap";
                    }
                    if (e[i].nbad != 0 && e[i].solvername != "mptiluc") {
                        str += (eng) ? " (or choose more reliable linear solver, for example, MPTILUC)" :
                                       " (или выбрать более надежный решатель, например, MPTILUC)";
                    }
                    str += ".]\n";
                }
            }
            if (str != "")
                str = (eng ? "Recommendation on linear solver parameteres:\n" :
                             "Рекомендации по параметрам линейного решателя:\n") + str;
            return str;
        }

#endif //USE_SOLVER_STAT

#if defined(USE_OPTIMIZER)
        struct TTSPOptions {
            bool enabled;
        };

        static TTSPConfiguration                                   *configuration = nullptr;
        static std::map<std::string, TTSPOptions>                  g_options      = std::map<std::string, TTSPOptions>();
        static std::map<std::string, INMOST::OptimizerInterface *> g_optimizers   = std::map<std::string, INMOST::OptimizerInterface *>();
#endif

        void Initialize(const std::string &path) {
			(void)path;
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
			(void)solver_name,(void)solver_prefix;
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
			(void)solver_name,(void)solver_prefix;
#if defined(USE_OPTIMIZER)
            const std::string &key = solver_name + ":" + solver_prefix;
            if (g_options.find(key) != g_options.end()) {
                g_options[key].enabled = false;
            }
#endif
        }

        bool isEnabled(const std::string &solver_name, const std::string &solver_prefix) {
			(void)solver_name,(void)solver_prefix;
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
			(void)solver_name,(void)solver_prefix;
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
			(void)solver,(void)use_last_suggestion;
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
			(void)solver,(void)metrics,(void)is_good;
#if defined(USE_SOLVER_STAT)
            SolverOptimizeSaveResult3(solver, metrics, is_good);
#endif
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
			(void)solver;
#if defined(USE_OPTIMIZER)
            const std::string &key = solver.SolverName() + ":" + solver.SolverPrefix();
            if (g_optimizers.find(key) != g_optimizers.end()) {
                delete g_optimizers[key];
                g_optimizers.erase(key);
            }
#endif
        }

        std::string SolutionMetadataLine() {
#if defined(USE_SOLVER_STAT)
            return SolutionMetadataLine3();
#else
            return "";
#endif
        }

    }

}
#endif // USE_SOLVER
