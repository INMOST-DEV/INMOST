//
// Created by Dmitri Bagaev on 2019-03-22.
//

#include "bayesian.h"

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
#include <limbo/model/multi_gp.hpp>
#include <limbo/acqui/gp_ucb.hpp>
#include <limbo/acqui/ei.hpp>

namespace INMOST {

    BayesianUniformDistribution::BayesianUniformDistribution() : distribution(0.0, 1.0) {
        int rank, size;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        unsigned int seed = static_cast<unsigned int>(time(NULL));

        if (rank == 0) {
            for (int i = 1; i < size; i++) {
                MPI_Send(&seed, 1, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD);
            }
        } else {
            MPI_Status status;
            MPI_Recv(&seed, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, &status);
        }

        generator.seed(seed);
    }

    double BayesianUniformDistribution::next() {
        return distribution(generator);
    }

    unsigned int BayesianUniformDistribution::operator()(unsigned int n) {
        return static_cast<unsigned int>(next() * n);
    }

    unsigned int    BayesianOptimizer::DEFAULT_UNIQUE_POINTS_MAX_COUNT    = 7;
    unsigned int    BayesianOptimizer::DEFAULT_UNIQUE_POINTS_RANDOM_COUNT = 5;
    unsigned int    BayesianOptimizer::DEFAULT_INITIAL_ITERATIONS_COUNT   = 5;
    double          BayesianOptimizer::DEFAULT_INITIAL_ITERATIONS_RADIUS  = 0.2;
    double          BayesianOptimizer::DEFAULT_MAX_JUMP_BARRIER           = 0.1;

    BayesianOptimizer::BayesianOptimizer(const std::string &name, const OptimizationParameters &space, const OptimizerProperties &properties, std::size_t buffer_capacity) :
            OptimizerInterface(name, space, properties, buffer_capacity),
            unique_points_max_count(BayesianOptimizer::DEFAULT_UNIQUE_POINTS_MAX_COUNT),
            unique_points_random_count(BayesianOptimizer::DEFAULT_UNIQUE_POINTS_RANDOM_COUNT),
            initial_iterations_count(BayesianOptimizer::DEFAULT_INITIAL_ITERATIONS_COUNT),
            initial_iterations_radius(BayesianOptimizer::DEFAULT_INITIAL_ITERATIONS_RADIUS),
            max_jump_barrier(BayesianOptimizer::DEFAULT_MAX_JUMP_BARRIER) {

        if (this->HasProperty("unique_points_max_count")) {
            unique_points_max_count = static_cast<unsigned int>(std::atoi(this->GetProperty("unique_points_max_count").c_str()));
        }

        if (this->HasProperty("unique_points_random_count")) {
            unique_points_random_count = static_cast<unsigned int>(std::atoi(this->GetProperty("unique_points_random_count").c_str()));
        }

        if (this->HasProperty("initial_iterations_count")) {
            initial_iterations_count = static_cast<unsigned int>(std::atoi(this->GetProperty("initial_iterations_count").c_str()));
        }

        if (this->HasProperty("initial_iterations_radius")) {
            initial_iterations_radius = std::atof(this->GetProperty("initial_iterations_radius").c_str());
        }

        if (this->HasProperty("max_jump_barrier")) {
            max_jump_barrier = std::atof(this->GetProperty("max_jump_barrier").c_str());
        }

    }

    SuggestionChangedParameters BayesianOptimizer::AlgorithmMakeSuggestion(const std::function<OptimizationFunctionInvokeResult(const OptimizationParameterPoints &,
                                                                                                                                const OptimizationParameterPoints &,
                                                                                                                                void *)> &invoke, void *data) const {

        auto unique  = results.GetLastUniqueEntries(unique_points_max_count);
        auto entries = parameters.GetParameterEntries();

        std::vector<SuggestionChangedParameter> changed_parameters;

        if (unique.size() == 0) {

            std::size_t index = 0;
            std::for_each(entries.cbegin(), entries.cend(), [&index, &changed_parameters](const OptimizationParametersEntry &entry) {
                changed_parameters.emplace_back(SuggestionChangedParameter(index++, entry.first.GetDefaultValue()));
            });

            return changed_parameters;
        } else if (unique.size() < initial_iterations_count) {

            std::size_t index = 0;
            std::for_each(entries.cbegin(), entries.cend(), [&index, &changed_parameters, this, &unique](const OptimizationParametersEntry &entry) {
                auto parameter = entry.first;

                double min_bound = parameter.GetMinimalValue();
                double max_bound = parameter.GetMaximumValue();

                double r = random.next() * (max_bound - min_bound) * initial_iterations_radius;

                double next = entry.second + r * (2.0 * (unique.size() % 2) - 1.0);

                while (next < min_bound || next > max_bound) {
                    r    = random.next() * (max_bound - min_bound) * initial_iterations_radius;
                    next = entry.second + r * (2.0 * (unique.size() % 2) - 1.0);
                }

                changed_parameters.emplace_back(SuggestionChangedParameter(index++, next));
            });

            return changed_parameters;
        }

        std::random_shuffle(unique.begin() + 1, unique.end(), random);

        struct Params {
            struct kernel {
                // BO_PARAM(double, noise, 0.01); default
                // BO_PARAM(bool, optimize_noise, false); default
                BO_PARAM(double, noise, 0.01);

                BO_PARAM(bool, optimize_noise, false);
            };

            struct kernel_squared_exp_ard : public limbo::defaults::kernel_squared_exp_ard {
                // BO_PARAM(int, k, 0); default
                // BO_PARAM(double, sigma_sq, 1); default

                BO_PARAM(int, k, 4);

                BO_PARAM(double, sigma_sq, 0.2);
            };

            struct opt_rprop : public limbo::defaults::opt_rprop {
            };
            struct opt_nloptnograd : public limbo::defaults::opt_nloptnograd {
            };

            struct acqui_ucb {
                // BO_PARAM(double, alpha, 0.5); default

                BO_PARAM(double, alpha, 0.25);
            };

            struct acqui_gpucb : public limbo::defaults::acqui_gpucb {
            };
            struct acqui_ei : public limbo::defaults::acqui_ei {
            };
        };

        using Kernel2_t = limbo::kernel::SquaredExpARD<Params>;
        using Mean_t = limbo::mean::FunctionARD<Params, limbo::mean::Data<Params>>;
        using GP2_t = limbo::model::GP<Params, Kernel2_t, Mean_t, limbo::model::gp::KernelMeanLFOpt<Params>>;

        GP2_t gp_ard;

        std::vector<Eigen::VectorXd> samples;
        std::vector<Eigen::VectorXd> observations;

        std::for_each(unique.cbegin(), unique.cbegin() + std::min(static_cast<std::size_t>(unique_points_random_count), unique.size()), [this, &samples, &observations]
                (const OptimizationParameterResult &result) {
            Eigen::VectorXd sample(parameters.Size());

            int i = 0;
            std::for_each(result.GetPointsAfter().cbegin(), result.GetPointsAfter().cend(), [&i, &sample, this](const OptimizationParameterPoint &point) {
                auto parameter = parameters.GetParameter(static_cast<size_t>(i));

                double min_bound = parameter.GetMinimalValue();
                double max_bound = parameter.GetMaximumValue();

                double normalized = (parameter.ExtractValueFromPoint(point) - min_bound) / (max_bound - min_bound); // Normalize here to [0, 1]

                sample(i) = normalized;

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
            return limbo::tools::make_vector((ob(0) - min_observation) / ((10.0 / 6.0) * (max_observation - min_observation)) + 0.2); // Normalize here to [ 0.2, 0.8 ]
        });

        gp_ard.compute(samples, observations);
        gp_ard.optimize_hyperparams();

        using acquiopt_t = limbo::opt::NLOptNoGrad<Params, nlopt::GN_ISRES>;
        using acquisition_function_t = limbo::acqui::UCB<Params, GP2_t>;

        acquiopt_t             acquiopt;
        acquisition_function_t acqui(gp_ard, 0);

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

            double barrier = max_jump_barrier;
            if (std::abs(new_sample(k) - starting_point(k)) > barrier) {
                new_sample(k) = starting_point(k) + barrier * (new_sample(k) - starting_point(k));
            }

            changed_parameters.emplace_back(SuggestionChangedParameter(k, (new_sample(k) * (max_bound - min_bound)) + min_bound));
        }

        return changed_parameters;
    }

    bool BayesianOptimizer::UpdateSpaceWithLatestResults() {
        const OptimizationParameterResult &last = results.at(0);


        if (last.IsGood() && (last.GetMetricsBefore() < 0.0 || (last.GetMetricsAfter() < last.GetMetricsBefore()))) {
            parameters.Update(last.GetChangedParameters(), last.GetMetricsAfter());
            return true;
        } else {
            return false;
        }
    }

    BayesianOptimizer::~BayesianOptimizer() {}


}
