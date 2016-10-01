#include "SolverDDPQILUC2.h"

namespace INMOST {

    SolverDDPQILUC2::SolverDDPQILUC2() {
        Method *preconditioner = new ILUC_preconditioner(info);
        solver = new KSOLVER(preconditioner, info);
        matrix = NULL;
    }

    SolverDDPQILUC2::SolverDDPQILUC2(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverDDPQILUC2::Assign(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverDDPQILUC2::Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix) {

    }

    void SolverDDPQILUC2::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) {
        if (matrix != NULL) {
            delete matrix;
        }
        matrix = new Sparse::Matrix(A);
        info.PrepareMatrix(*matrix, parameters.GetParameterEnum("additive_schwartz_overlap"));
        solver->ReplaceMAT(*matrix);

        solver->RealParameter(":tau") = parameters.GetParameterReal("drop_tolerance");
        solver->RealParameter(":tau2") = parameters.GetParameterReal("reuse_tolerance");
        solver->EnumParameter(":scale_iters") = parameters.GetParameterEnum("rescale_iterations");
        solver->RealParameter(":ddpq_tau") = parameters.GetParameterReal("ddpq_tolerance");
        solver->EnumParameter(":reorder_nnz") = parameters.GetParameterEnum("reorder_nonzero");
        solver->EnumParameter(":estimator") = parameters.GetParameterEnum("condition_estimation");
        solver->EnumParameter(":ddpq_tau_adapt") = parameters.GetParameterEnum("adapt_ddpq_tolerance");

        if (sizeof(KSOLVER) == sizeof(BCGSL_solver)) {
            solver->EnumParameter("levels") = parameters.GetParameterEnum("gmres_substeps");
        }

        if (!solver->isInitialized()) {
            solver->Initialize();
        }
    }

    bool SolverDDPQILUC2::Solve(Sparse::Vector &RHS, Sparse::Vector &SOL) {
        solver->EnumParameter("maxits") = parameters.GetParameterEnum("maximum_iterations");
        solver->RealParameter("rtol") = parameters.GetParameterReal("relative_tolerance");
        solver->RealParameter("atol") = parameters.GetParameterReal("absolute_tolerance");
        solver->RealParameter("divtol") = parameters.GetParameterReal("divergence_tolerance");

        return solver->Solve(RHS, SOL);
    }

    bool SolverDDPQILUC2::Clear() {
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

    bool SolverDDPQILUC2::isMatrixSet() {
        return matrix != NULL;
    }

    const INMOST_DATA_ENUM_TYPE SolverDDPQILUC2::Iterations() const {
        return solver->GetIterations();
    }

    const INMOST_DATA_REAL_TYPE SolverDDPQILUC2::Residual() const {
        return solver->GetResidual();
    }

    const std::string SolverDDPQILUC2::ReturnReason() const {
        return solver->GetReason();
    }

    const std::string SolverDDPQILUC2::SolverName() const {
        return "inner_ddpqiluc2";
    }

    void SolverDDPQILUC2::Finalize() {

    }

    SolverDDPQILUC2::~SolverDDPQILUC2() {
        this->Clear();
    }

}
