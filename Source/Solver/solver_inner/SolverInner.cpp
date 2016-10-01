#include "SolverInner.h"

namespace INMOST {

    SolverInner::SolverInner() {

    }

    SolverInner::SolverInner(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverInner::Assign(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverInner::Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix) {

    }

    bool SolverInner::Solve(Sparse::Vector &RHS, Sparse::Vector &SOL) {
        solver->EnumParameter("maxits") = parameters.GetParameterEnum("maximum_iterations");
        solver->RealParameter("rtol") = parameters.GetParameterReal("relative_tolerance");
        solver->RealParameter("atol") = parameters.GetParameterReal("absolute_tolerance");
        solver->RealParameter("divtol") = parameters.GetParameterReal("divergence_tolerance");

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

    void SolverInner::Finalize() {

    }

    SolverInner::~SolverInner() {
        this->Clear();
    }
}