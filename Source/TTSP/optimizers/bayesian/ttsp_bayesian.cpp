//
// Created by Dmitri Bagaev on 2019-03-22.
//

#include "ttsp_bayesian.h"

#define USE_NLOPT

#include <nlopt.hpp>

#include <boost/parameter/aux_/void.hpp>

#include <limbo/opt/nlopt_no_grad.hpp>
#include <Eigen/Core>
#include <limbo/kernel/exp.hpp>
#include <limbo/kernel/squared_exp_ard.hpp>
#include <limbo/mean/function_ard.hpp>
#include <limbo/model/gp.hpp>
#include <limbo/model/gp/kernel_lf_opt.hpp>
#include <limbo/tools.hpp>
#include <limbo/tools/macros.hpp>
#include <limbo/bayes_opt/bo_base.hpp>
#include <limbo/bayes_opt/boptimizer.hpp>
#include <limbo/model/gp/kernel_mean_lf_opt.hpp>
#include <limbo/acqui/gp_ucb.hpp>
#include <limbo/acqui/ei.hpp>

namespace TTSP {


    BayesianOptimizer::BayesianOptimizer(const OptimizationParameters &space, const OptimizerProperties &properties, std::size_t buffer_capacity) :
            OptimizerInterface(space, properties, buffer_capacity), current_iteration(0) {}

    OptimizationAlgorithmSuggestion BayesianOptimizer::AlgorithmMakeSuggestion(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                                                    const OptimizationParameterPoints &,
                                                                                                                                    void *)> &invoke, void *data) const {

        if (results.size() == 0) {
            return std::make_pair(0, parameters.GetParameter(0).GetDefaultValue());
        }

        struct Params {
            struct kernel : public limbo::defaults::kernel {
            };
            struct kernel_squared_exp_ard : public limbo::defaults::kernel_squared_exp_ard {
            };
            struct opt_rprop : public limbo::defaults::opt_rprop {
            };
            struct opt_nloptnograd : public limbo::defaults::opt_nloptnograd {
            };
            struct acqui_ucb : public limbo::defaults::acqui_ucb {
            };
            struct acqui_ei : public limbo::defaults::acqui_ei {
            };
        };

        using Kernel2_t = limbo::kernel::SquaredExpARD<Params>;
        using Mean_t = limbo::mean::FunctionARD<Params, limbo::mean::Data<Params>>;
        using GP2_t = limbo::model::GP<Params, Kernel2_t, Mean_t, limbo::model::gp::KernelMeanLFOpt<Params>>;

        GP2_t gp_ard;
        // do not forget to call the optimization!

        std::vector<Eigen::VectorXd> samples;
        std::vector<Eigen::VectorXd> observations;

        std::for_each(results.cbegin(), results.cend(), [this, &samples, &observations](const OptimizationParameterResult &result) {
            Eigen::VectorXd sample(parameters.Size());

            int i = 0;
            std::for_each(result.GetPointsAfter().cbegin(), result.GetPointsAfter().cend(), [&i, &sample, this](const OptimizationParameterPoint &point) {
                auto parameter = parameters.GetParameter(static_cast<size_t>(i));

                double min_bound = parameter.GetMinimalValue();
                double max_bound = parameter.GetMaximumValue();

                sample(i) = (parameter.ExtractValueFromPoint(point) - min_bound) / (max_bound - min_bound); // Normalize here to [0, 1]
                i += 1;
            });

            samples.push_back(sample);
            observations.push_back(limbo::tools::make_vector(-1.0 * result.GetMetricsAfter()));
        });

        double min_observation = (*std::min_element(observations.cbegin(), observations.cend(), [](const Eigen::VectorXd &l, const Eigen::VectorXd &r) {
            return l(0) < r(0);
        }))(0);

        double max_observation = (*std::max_element(observations.cbegin(), observations.cend(), [](const Eigen::VectorXd &l, const Eigen::VectorXd &r) {
            return l(0) < r(0);
        }))(0);

        std::transform(observations.cbegin(), observations.cend(), observations.begin(), [min_observation, max_observation](const Eigen::VectorXd &ob) {
            return limbo::tools::make_vector((ob(0) - min_observation) / (max_observation - min_observation));
        });


        gp_ard.compute(samples, observations);
        gp_ard.optimize_hyperparams();

        using acquiopt_t = limbo::opt::NLOptNoGrad<Params, nlopt::GN_DIRECT_L_RAND>;
        using acquisition_function_t = limbo::acqui::UCB<Params, GP2_t>;

        acquiopt_t             acquiopt;
        acquisition_function_t acqui(gp_ard, static_cast<int>(current_iteration));

        auto afun               = limbo::FirstElem();
        auto acqui_optimization = [&](const Eigen::VectorXd &x, bool g) { return acqui(x, afun, g); };

        Eigen::VectorXd starting_point(parameters.Size());

        int i = 0;
        std::for_each(parameters.GetParameterEntries().cbegin(), parameters.GetParameterEntries().cend(), [&i, &starting_point](const OptimizationParametersEntry &entry) {

            auto parameter = entry.first;

            double min_bound = parameter.GetMinimalValue();
            double max_bound = parameter.GetMaximumValue();

            starting_point(i) = (entry.second - min_bound) / (max_bound - min_bound);

            i += 1;
        });

        Eigen::VectorXd new_sample = acquiopt(acqui_optimization, starting_point, true);

        for (int k = 0; k < parameters.Size(); ++k) {
            auto parameter = parameters.GetParameter(static_cast<size_t>(k));

            double min_bound = parameter.GetMinimalValue();
            double max_bound = parameter.GetMaximumValue();

            new_sample(k) = (new_sample(k) * (max_bound - min_bound)) + min_bound;
        }

        current_iteration += 1;

        return std::make_pair(0, new_sample[0]);
    }

    BayesianOptimizer::~BayesianOptimizer() {}


}
