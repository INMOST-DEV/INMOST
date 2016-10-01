#include "SolverMPTILUC.h"

namespace INMOST {

    SolverMPTILUC::SolverMPTILUC(){
        Method *preconditioner = new MTILUC_preconditioner(info);
        solver = new KSOLVER(preconditioner, info);
        matrix = NULL;
    }

    SolverMPTILUC::SolverMPTILUC(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverMPTILUC::Assign(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverMPTILUC::Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix) {

    }

    void SolverMPTILUC::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) {
        if (matrix != NULL) {
            delete matrix;
        }
        matrix = new Sparse::Matrix(A);
        info.PrepareMatrix(*matrix, parameters.GetParameterEnum("additive_schwartz_overlap"));
        solver->ReplaceMAT(*matrix);

        solver->RealParameter(":tau") = parameters.GetParameterReal("drop_tolerance");
        solver->RealParameter(":tau2") = parameters.GetParameterReal("reuse_tolerance");
        solver->EnumParameter(":scale_iters") = parameters.GetParameterEnum("rescale_iterations");
        solver->EnumParameter(":estimator") = parameters.GetParameterEnum("condition_estimation");

        if (sizeof(KSOLVER) == sizeof(BCGSL_solver)) {
            solver->EnumParameter("levels") = parameters.GetParameterEnum("gmres_substeps");
        }

        if (!solver->isInitialized()) {
            solver->Initialize();
        }
    }

    bool SolverMPTILUC::Solve(Sparse::Vector &RHS, Sparse::Vector &SOL) {
        solver->EnumParameter("maxits") = parameters.GetParameterEnum("maximum_iterations");
        solver->RealParameter("rtol") = parameters.GetParameterReal("relative_tolerance");
        solver->RealParameter("atol") = parameters.GetParameterReal("absolute_tolerance");
        solver->RealParameter("divtol") = parameters.GetParameterReal("divergence_tolerance");

        return solver->Solve(RHS, SOL);
    }

    bool SolverMPTILUC::Clear() {
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

    bool SolverMPTILUC::isMatrixSet() {
        return matrix != NULL;
    }

    const INMOST_DATA_ENUM_TYPE SolverMPTILUC::Iterations() const {
        return solver->GetIterations();
    }

    const INMOST_DATA_REAL_TYPE SolverMPTILUC::Residual() const {
        return solver->GetResidual();
    }

    const std::string SolverMPTILUC::ReturnReason() const {
        return solver->GetReason();
    }

    const std::string SolverMPTILUC::SolverName() const {
        return "inner_mptiluc";
    }

    void SolverMPTILUC::Finalize() {

    }

    SolverMPTILUC::~SolverMPTILUC() {
        this->Clear();
    }

}
