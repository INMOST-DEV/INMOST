#include "SolverInner.h"

namespace INMOST {

    SolverInner::SolverInner(SolverParameters &parameters): SolverInterface(parameters) {
        parameters.require("additive_schwartz_overlap", "1");
        parameters.require("maximum_iterations", "2500");
        parameters.require("reorder_nonzero", "1");
        parameters.require("rescale_iterations", "6");
        parameters.require("condition_estimation", "1");
        parameters.require("adapt_ddpq_tolerance", "1");
        parameters.require("gmres_substeps", "2");

        parameters.require("absolute_tolerance", "1.0e-5");
        parameters.require("relative_tolerance", "1.0e-12");
        parameters.require("divergence_tolerance", "1.0e+100");
        parameters.require("drop_tolerance", "0.005");
        parameters.require("reuse_tolerance", "0.00005");
        parameters.require("ddpq_tolerance", "0.75");
        parameters.require("fill_level", "3");
    }

    SolverInner::SolverInner(const SolverInterface *other): SolverInterface(other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverInner::Assign(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverInner::Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix) {
//        FILE *databaseFile = fopen(parameters_file, "r");
//        if (!databaseFile) {
//            return;
//        }
//        char *tmp = (char *) calloc(256, sizeof(char));
//        char *parameterName = (char *) calloc(128, sizeof(char));
//        char *parameterValue = (char *) calloc(128, sizeof(char));
//        while (!feof(databaseFile) && fgets(tmp, 256, databaseFile)) {
//            char *line = tmp;
//            //Comment str
//            if (line[0] == '#') continue;
//            sscanf(line, "%s %s", parameterName, parameterValue);
//            parameters.require(parameterName, parameterValue);
//        }
//        free(parameterValue);
//        free(parameterName);
//        free(tmp);
    }

    bool SolverInner::Solve(Sparse::Vector &RHS, Sparse::Vector &SOL) {
        solver->EnumParameter("maxits") = parameters.get < INMOST_DATA_ENUM_TYPE > ("maximum_iterations");
        solver->RealParameter("rtol") = parameters.get<INMOST_DATA_REAL_TYPE>("relative_tolerance");
        solver->RealParameter("atol") = parameters.get<INMOST_DATA_REAL_TYPE>("absolute_tolerance");
        solver->RealParameter("divtol") = parameters.get<INMOST_DATA_REAL_TYPE>("divergence_tolerance");

        return solver->Solve(RHS, SOL);
    }

    bool SolverInner::Clear() {
        info.Clear();
        if (matrix != NULL) {
            delete matrix;
            matrix = NULL;
        }
        if (solver != NULL) {
            delete solver;
            solver = NULL;
        }
        return true;
    }

    bool SolverInner::isMatrixSet() {
        return matrix != NULL;
    }

    const INMOST_DATA_ENUM_TYPE SolverInner::Iterations() const {
        return solver->GetIterations();
    }

    const INMOST_DATA_REAL_TYPE SolverInner::Residual() const {
        return solver->GetResidual();
    }

    const std::string SolverInner::ReturnReason() const {
        return solver->GetReason();
    }

    const INMOST_DATA_REAL_TYPE SolverInner::Condest(INMOST_DATA_REAL_TYPE tol, INMOST_DATA_ENUM_TYPE maxiter) {
#if defined(ACCELERATED_CONDEST)
        INMOST_DATA_ENUM_TYPE lbeg, lend, l, iter;
        INMOST_DATA_REAL_TYPE norm, sum[2], norm_prev, lambda_min, lambda_max;
        bool diverged_max = false, diverged_min = false;
        info.GetLocalRegion(info.GetRank(), lbeg, lend);
        Sparse::Vector v, Av;
        info.PrepareVector(v);
        info.PrepareVector(Av);
        //Set v to random
        norm = 0;
        for (l = lbeg; l < lend; ++l) {
            v[l] = rand() / static_cast<INMOST_DATA_REAL_TYPE>(RAND_MAX);
            norm += v[l] * v[l];
        }
        info.Integrate(&norm, 1);
        norm_prev = 0.0;
        norm = sqrt(norm);
        for (l = lbeg; l < lend; ++l) v[l] /= norm;
        //Compute maximum eigenvalue
        iter = 0;
        while (fabs((norm - norm_prev) / norm) > tol && iter < maxiter) {
            info.Update(v);
#if defined(USE_OMP)
#pragma omp parallel
#endif
            matrix->MatVec(1.0, v, 0.0, Av);
            norm_prev = norm;
            sum[0] = sum[1] = 0.0;
            for (l = lbeg; l < lend; ++l) {
                sum[0] += v[l] * Av[l];
                sum[1] += v[l] * v[l];
            }
            info.Integrate(sum, 2);
            norm = fabs(sum[0]) / sum[1];
            for (l = lbeg; l < lend; ++l) v[l] = Av[l] / norm;
#if defined(PRINT_CONDEST)
            std::cout << "iteration " << iter << " norm " << norm << std::endl;
#endif
            iter++;
        }
#if defined(PRINT_CONDEST)
        std::cout << "lambda_max " << norm << std::endl;
#endif
        if (iter == maxiter) {
            diverged_max = true;
            std::cout << "Max not converged" << std::endl;
        }
        lambda_max = norm;
        //Set v to random
        norm = 0;
        for (l = lbeg; l < lend; ++l) {
            v[l] = rand() / static_cast<INMOST_DATA_REAL_TYPE>(RAND_MAX);
            norm += v[l] * v[l];
        }
        info.Integrate(&norm, 1);
        norm_prev = 0.0;
        norm = sqrt(norm);
        for (l = lbeg; l < lend; ++l) v[l] /= norm;
        //Compute minimal eigenvalue
        iter = 0;
        while (fabs((norm - norm_prev) / norm) > tol && iter < maxiter) {
            info.Update(v);
            Solve(v, Av);
            norm_prev = norm;
            sum[0] = sum[1] = 0;
            for (l = lbeg; l < lend; ++l) {
                sum[0] += v[l] * Av[l];
                sum[1] += v[l] * v[l];
            }
            info.Integrate(sum, 2);
            norm = fabs(sum[0]) / sum[1];
            for (l = lbeg; l < lend; ++l) v[l] = Av[l] / norm;
#if defined(PRINT_CONDEST)
            std::cout << "iteration " << iter << " norm " << norm << "\t\t" << std::endl;
#endif
            iter++;
        }
#if defined(PRINT_CONDEST)
        std::cout << "lambda_min " << 1.0 / norm << std::endl;
#endif
        if (iter == maxiter) {
            diverged_min = true;
            std::cout << "Min not converged" << std::endl;
        }
        lambda_min = 1.0 / norm;
        if (diverged_max || diverged_min)
            return 1.0e+100;
#if defined(PRINT_CONDEST)
        std::cout << "Condest: " << lambda_max / lambda_min << std::endl;
#endif
        return lambda_max / lambda_min;
#else //!ACCELERATED_CONDEST
        INMOST_DATA_ENUM_TYPE lbeg, lend, l, iter;
        INMOST_DATA_REAL_TYPE norm, norm_prev, lambda_min, lambda_max;
        bool diverged_max = false, diverged_min = false;
        info.GetLocalRegion(info.GetRank(), lbeg, lend);
        Sparse::Vector v, Av;
        info.PrepareVector(v);
        info.PrepareVector(Av);
        //Set v to random
        norm = 0;
        for (l = lbeg; l < lend; ++l) {
            v[l] = rand() / static_cast<INMOST_DATA_REAL_TYPE>(RAND_MAX);
            norm += v[l] * v[l];
        }
        info.Integrate(&norm, 1);
        norm_prev = 0.0;
        norm = sqrt(norm);
        for (l = lbeg; l < lend; ++l) v[l] /= norm;
        //Compute maximum eigenvalue
        iter = 0;
        while (fabs(norm - norm_prev) / norm > tol && iter < maxiter) {
            info.Update(v);
#if defined(USE_OMP)
#pragma omp parallel
#endif
            matrix->MatVec(1.0, v, 0.0, Av);
            v.Swap(Av);
            norm_prev = norm;
            norm = 0.0;
            for (l = lbeg; l < lend; ++l) norm += v[l] * v[l];
            info.Integrate(&norm, 1);
            norm = sqrt(norm);
            for (l = lbeg; l < lend; ++l) v[l] /= norm;
#if defined(PRINT_CONDEST)
            std::cout << "iteration " << iter << " norm " << norm << std::endl;
#endif
            iter++;
        }
#if defined(PRINT_CONDEST)
        std::cout << "lambda_max " << norm << std::endl;
#endif
        if (iter == maxiter) {
            norm = std::max(norm, norm_prev);
            //diverged_max = true;
            //std::cout << "Max not converged" << std::endl;
        }
        lambda_max = norm;
        //Set v to random
        norm = 0;
        for (l = lbeg; l < lend; ++l) {
            v[l] = rand() / static_cast<INMOST_DATA_REAL_TYPE>(RAND_MAX);
            norm += v[l] * v[l];
        }
        info.Integrate(&norm, 1);
        norm_prev = 0.0;
        norm = sqrt(norm);
        for (l = lbeg; l < lend; ++l) v[l] /= norm;
        //Compute minimal eigenvalue
        iter = 0;
        while (fabs(norm - norm_prev) / norm > tol && iter < maxiter) {
            info.Update(v);
            Solve(v, Av);
            v.Swap(Av);
            norm_prev = norm;
            norm = 0.0;
            for (l = lbeg; l < lend; ++l) norm += v[l] * v[l];
            info.Integrate(&norm, 1);
            norm = sqrt(norm);
            for (l = lbeg; l < lend; ++l) v[l] /= norm;
#if defined(PRINT_CONDEST)
            std::cout << "iteration " << iter << " norm " << norm << "\t\t" << std::endl;
#endif
            iter++;
        }
#if defined(PRINT_CONDEST)
        std::cout << "lambda_min " << 1.0 / norm << std::endl;
#endif
        if (iter == maxiter) {
            norm = std::max(norm, norm_prev);
            //diverged_min = true;
            //std::cout << "Min not converged" << std::endl;
        }
        lambda_min = 1.0 / norm;
        //Condition is always false? Maybe check this
        if (diverged_max || diverged_min)
            return 1.0e+100;
#if defined(PRINT_CONDEST)
        std::cout << "Condest: " << lambda_max / lambda_min << std::endl;
#endif
        return lambda_max / lambda_min;
#endif
    }

    void SolverInner::Finalize() {

    }

    SolverInner::~SolverInner() {
        this->Clear();
    }

}