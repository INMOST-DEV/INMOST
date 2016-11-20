#include "SolverInner.h"

namespace INMOST {

    SolverInner::SolverInner() {
        iters = 2500;
        scale_iters = 6;
        condest = 1;
        ddpqatol = 1;
        overlap = 1;
        ell = 2;
        reorder_nnz = 1;

        atol = 1.0e-5;
        rtol = 1.0e-12;
        dtol = 1.0e+100;
        tau = 0.005;
        tau2 = 0.00005;
        ddpqtol = 0.75;
        fill = 3;
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
        solver->EnumParameter("maxits") = iters;
        solver->RealParameter("rtol") = rtol;
        solver->RealParameter("atol") = atol;
        solver->RealParameter("divtol") = dtol;

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

    std::string SolverInner::GetParameter(std::string name) const {
        if(name == "maximum_iterations" ) return to_string(iters);
        else if( name == "rescale_iterations" ) return to_string(scale_iters);
        else if( name == "condition_estimation" ) return to_string(condest);
        else if( name == "adapt_ddpq_tolerance" ) return to_string(ddpqatol);
        else if( name == "schwartz_overlap" ) return to_string(overlap);
        else if( name == "gmres_substeps" ) return to_string(ell);
        else if( name == "reorder_nonzeros" ) return to_string(reorder_nnz);
        else if( name == "absolute_tolerance") return to_string(atol);
        else if( name == "relative_tolerance") return to_string(rtol);
        else if( name == "divergence_tolerance") return to_string(dtol);
        else if( name == "drop_tolerance") return to_string(tau);
        else if( name == "reuse_tolerance") return to_string(tau2);
        else if( name == "ddpq_tolerance") return to_string(ddpqtol);
        else if( name == "fill_level") return to_string(fill);
        else {
            std::cout << "Parameter " << name << " is unknown" << std::endl;
            return "";
        }
    }

    void SolverInner::SetParameter(std::string name, std::string value) {
        const char *val = value.c_str();
        if(name == "maximum_iterations" ) iters = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if( name == "rescale_iterations" ) scale_iters = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if( name == "condition_estimation" ) condest = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if( name == "adapt_ddpq_tolerance" ) ddpqatol = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if( name == "schwartz_overlap" ) overlap = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if( name == "gmres_substeps" ) ell = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if( name == "reorder_nonzeros" ) reorder_nnz = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if( name == "absolute_tolerance") atol = atof(val);
        else if( name == "relative_tolerance") rtol = atof(val);
        else if( name == "divergence_tolerance") dtol = atof(val);
        else if( name == "drop_tolerance") tau = atof(val);
        else if( name == "reuse_tolerance") tau2  = atof(val);
        else if( name == "ddpq_tolerance") ddpqtol = atof(val);
        else if( name == "fill_level") fill = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else std::cout << "Parameter " << name << " is unknown" << std::endl;
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