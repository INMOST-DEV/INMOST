#include "SolverDDPQILUC2.h"

namespace INMOST {

    SolverDDPQILUC2::SolverDDPQILUC2() {
        Method *preconditioner = new ILUC_preconditioner(info);
        solver = new KSOLVER(preconditioner, info);
        matrix = NULL;
        additive_schwartz_overlap = DEFAULT_ADDITIVE_SCHWARTZ_OVERLAP;
        maximum_iterations = DEFAULT_MAXIMUM_ITERATIONS;
        absolute_tolerance = DEFAULT_ABSOLUTE_TOLERANCE;
        relative_tolerance = DEFAULT_RELATIVE_TOLERANCE;
        divergence_tolerance = DEFAULT_DIVERGENCE_TOLERANCE;

        preconditioner_drop_tolerance = DEFAULT_PRECONDITIONER_DROP_TOLERANCE;
        preconditioner_reuse_tolerance = DEFAULT_PRECONDITIONER_REUSE_TOLERANCE;
        preconditioner_rescale_iterations = DEFAULT_PRECONDITIONER_RESCALE_ITERS;
        preconditioner_ddpq_tolerance = DEFAULT_PRECONDITIONER_DDPQ_TOLERANCE;
        preconditioner_reorder_nonzero = DEFAULT_PRECONDITIONER_REORDER_NONZEROS;
        preconditioner_condition_estimation = DEFAULT_PRECONDITIONER_CONDITION_ESTIMATION;
        preconditioner_adapt_ddpq_tolerance = DEFAULT_PRECONDITIONER_ADAPT_DDPQ_TOLERANCE;
        solver_gmres_substeps = DEFAULT_SOLVER_GMRES_SUBSTEPS;
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
        info.PrepareMatrix(*matrix, additive_schwartz_overlap);
        solver->ReplaceMAT(*matrix);

        solver->RealParameter(":tau") = preconditioner_drop_tolerance;
        solver->RealParameter(":tau2") = preconditioner_reuse_tolerance;
        solver->EnumParameter(":scale_iters") = preconditioner_rescale_iterations;
        solver->RealParameter(":ddpq_tau") = preconditioner_ddpq_tolerance;
        solver->EnumParameter(":reorder_nnz") = preconditioner_reorder_nonzero;
        solver->EnumParameter(":estimator") = preconditioner_condition_estimation;
        solver->EnumParameter(":ddpq_tau_adapt") = preconditioner_adapt_ddpq_tolerance;

        if (sizeof(KSOLVER) == sizeof(BCGSL_solver)) {
            solver->EnumParameter("levels") = solver_gmres_substeps;
        }

        if (!solver->isInitialized()) {
            solver->Initialize();
        }
    }

    bool SolverDDPQILUC2::Solve(Sparse::Vector &RHS, Sparse::Vector &SOL) {
        solver->EnumParameter("maxits") = maximum_iterations;
        solver->RealParameter("rtol") = relative_tolerance;
        solver->RealParameter("atol") = absolute_tolerance;
        solver->RealParameter("divtol") = divergence_tolerance;

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

    INMOST_DATA_REAL_TYPE SolverDDPQILUC2::GetPropertyReal(std::string property) const {
        return solver->RealParameter(property);
    }

    INMOST_DATA_ENUM_TYPE SolverDDPQILUC2::GetPropertyEnum(std::string property) const {
        return solver->EnumParameter(property);
    }

    void SolverDDPQILUC2::SetPropertyReal(std::string property, INMOST_DATA_REAL_TYPE value) {
        if (property == "absolute_tolerance")
            absolute_tolerance = value;
        else if (property == "relative_tolerance")
            relative_tolerance = value;
        else if (property == "divergence_tolerance")
            divergence_tolerance = value;
        else if (property == "drop_tolerance")
            preconditioner_drop_tolerance = value;
        else if (property == "reuse_tolerance")
            preconditioner_reuse_tolerance = value;
        else if (property == "ddpq_tolerance")
            preconditioner_ddpq_tolerance = value;
        else std::cout << "Parameter " << property << " of real type is unknown" << std::endl;
    }

    void SolverDDPQILUC2::SetPropertyEnum(std::string property, INMOST_DATA_ENUM_TYPE value) {
        if (property == "maximum_iterations")
            maximum_iterations = value;
        else if (property == "rescale_iterations")
            preconditioner_rescale_iterations = value;
        else if (property == "condition_estimation")
            preconditioner_condition_estimation = value;
        else if (property == "adapt_ddpq_tolerance")
            preconditioner_adapt_ddpq_tolerance = value;
        else if (property == "schwartz_overlap")
            additive_schwartz_overlap = value;
        else if (property == "gmres_substeps")
            solver_gmres_substeps = value;
        else if (property == "reorder_nonzeros")
            preconditioner_reorder_nonzero = value;
        else std::cout << "Parameter " << property << " of integral type is unknown" << std::endl;
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
