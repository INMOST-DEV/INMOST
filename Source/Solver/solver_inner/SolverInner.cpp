#define _CRT_SECURE_NO_WARNINGS
#include "SolverInner.h"

namespace INMOST {

    SolverInner::SolverInner() {
        maximum_iterations = 2500;
        atol = 1.0e-5;
        rtol = 1.0e-12;
        dtol = 1.0e+100;
    }

    SolverInterface *SolverInner::Copy(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
        (void) other;
    }

    void SolverInner::Assign(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
        (void) other;
    }

    void SolverInner::Setup(int *argc, char ***argv, SolverParameters &p)
	{
        if (!p.internalFile.empty()) {
            SolverParameters::SetInnerParametersFromFile(p.internalFile, this);
        } else {
            for (parameters_iterator_t parameter = p.parameters.begin(); parameter < p.parameters.end(); parameter++)
                this->SetParameter((*parameter).first, (*parameter).second);
        }
        (void) argc; (void) argv;
    }

    bool SolverInner::Solve(Sparse::Vector &RHS, Sparse::Vector &SOL) {
        solver->EnumParameter("maxits") = maximum_iterations;
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
        if (name == "maximum_iterations") return to_string(maximum_iterations);
        else if (name == "absolute_tolerance") return to_string(atol);
        else if (name == "relative_tolerance") return to_string(rtol);
        else if (name == "divergence_tolerance") return to_string(dtol);
        else {
#if !defined(SILENCE_SET_PARAMETER)
            std::cout << "Parameter " << name << " is unknown" << std::endl;
#endif
            return "";
        }
    }

    void SolverInner::SetParameter(std::string name, std::string value) {
        const char *val = value.c_str();
        if (name == "maximum_iterations") maximum_iterations = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "absolute_tolerance") atol = atof(val);
        else if (name == "relative_tolerance") rtol = atof(val);
        else if (name == "divergence_tolerance") dtol = atof(val);
#if !defined(SILENCE_SET_PARAMETER)
        else std::cout << "Parameter " << name << " is unknown" << std::endl;
#endif
    }

    INMOST_DATA_ENUM_TYPE SolverInner::Iterations() const {
        return solver->GetIterations();
    }

    INMOST_DATA_REAL_TYPE SolverInner::Residual() const {
        return solver->GetResidual();
    }

    const std::string SolverInner::ReturnReason() const {
        return solver->GetReason();
    }

    INMOST_DATA_REAL_TYPE SolverInner::Condest(INMOST_DATA_REAL_TYPE tol, INMOST_DATA_ENUM_TYPE maxiter) {
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
