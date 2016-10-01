#include "SolverMPTILU2.h"

namespace INMOST {

    SolverMPTILU2::SolverMPTILU2() {
        Method *preconditioner = new MTILU2_preconditioner(info);
        solver = new KSOLVER(preconditioner, info);
        matrix = NULL;
    }

    SolverMPTILU2::SolverMPTILU2(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverMPTILU2::Assign(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverMPTILU2::Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix) {

    }

    void SolverMPTILU2::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) {
        if (matrix != NULL) {
            delete matrix;
        }
        matrix = new Sparse::Matrix(A);
        info.PrepareMatrix(*matrix, parameters.GetParameterEnum("additive_schwartz_overlap"));
        solver->ReplaceMAT(*matrix);

        solver->RealParameter(":tau") = parameters.GetParameterReal("drop_tolerance");
        solver->RealParameter(":tau2") = parameters.GetParameterReal("reuse_tolerance");
        solver->EnumParameter(":scale_iters") = parameters.GetParameterEnum("rescale_iterations");

        if (sizeof(KSOLVER) == sizeof(BCGSL_solver)) {
            solver->EnumParameter("levels") = parameters.GetParameterEnum("gmres_substeps");
        }

        if (!solver->isInitialized()) {
            solver->Initialize();
        }
    }

    bool SolverMPTILU2::Solve(Sparse::Vector &RHS, Sparse::Vector &SOL) {
        solver->EnumParameter("maxits") = parameters.GetParameterEnum("maximum_iterations");
        solver->RealParameter("rtol") = parameters.GetParameterReal("relative_tolerance");
        solver->RealParameter("atol") = parameters.GetParameterReal("absolute_tolerance");
        solver->RealParameter("divtol") = parameters.GetParameterReal("divergence_tolerance");

        return solver->Solve(RHS, SOL);
    }

    bool SolverMPTILU2::Clear() {
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

    bool SolverMPTILU2::isMatrixSet() {
        return matrix != NULL;
    }

    const INMOST_DATA_ENUM_TYPE SolverMPTILU2::Iterations() const {
        return solver->GetIterations();
    }

    const INMOST_DATA_REAL_TYPE SolverMPTILU2::Residual() const {
        return solver->GetResidual();
    }

    const std::string SolverMPTILU2::ReturnReason() const {
        return solver->GetReason();
    }

    const std::string SolverMPTILU2::SolverName() const {
        return "inner_mptilu2";
    }

    void SolverMPTILU2::Finalize() {

    }

    SolverMPTILU2::~SolverMPTILU2() {
        this->Clear();
    }

}
