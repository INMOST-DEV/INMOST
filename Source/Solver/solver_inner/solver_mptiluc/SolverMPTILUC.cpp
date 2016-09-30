#include "SolverMPTILUC.h"

namespace INMOST {

    SolverMPTILUC::SolverMPTILUC(){
        Method *preconditioner = new MTILUC_preconditioner(info);
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
        preconditioner_condition_estimation = DEFAULT_PRECONDITIONER_CONDITION_ESTIMATION;
        solver_gmres_substeps = DEFAULT_SOLVER_GMRES_SUBSTEPS;
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
        info.PrepareMatrix(*matrix, additive_schwartz_overlap);
        solver->ReplaceMAT(*matrix);

        solver->RealParameter(":tau") = preconditioner_drop_tolerance;
        solver->RealParameter(":tau2") = preconditioner_reuse_tolerance;
        solver->EnumParameter(":scale_iters") = preconditioner_rescale_iterations;
        solver->EnumParameter(":estimator") = preconditioner_condition_estimation;

        if (sizeof(KSOLVER) == sizeof(BCGSL_solver)) {
            solver->EnumParameter("levels") = solver_gmres_substeps;
        }

        if (!solver->isInitialized()) {
            solver->Initialize();
        }
    }

    bool SolverMPTILUC::Solve(Sparse::Vector &RHS, Sparse::Vector &SOL) {
        solver->EnumParameter("maxits") = maximum_iterations;
        solver->RealParameter("rtol") = relative_tolerance;
        solver->RealParameter("atol") = absolute_tolerance;
        solver->RealParameter("divtol") = divergence_tolerance;

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

    INMOST_DATA_REAL_TYPE SolverMPTILUC::GetPropertyReal(std::string property) const {
        return solver->RealParameter(property);
    }

    INMOST_DATA_ENUM_TYPE SolverMPTILUC::GetPropertyEnum(std::string property) const {
        return solver->EnumParameter(property);
    }

    void SolverMPTILUC::SetPropertyReal(std::string property, INMOST_DATA_REAL_TYPE value) {
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
        else std::cout << "Parameter " << property << " of real type is unknown" << std::endl;
    }

    void SolverMPTILUC::SetPropertyEnum(std::string property, INMOST_DATA_ENUM_TYPE value) {
        if (property == "maximum_iterations")
            maximum_iterations = value;
        else if (property == "rescale_iterations")
            preconditioner_rescale_iterations = value;
        else if (property == "condition_estimation")
            preconditioner_condition_estimation = value;
        else if (property == "schwartz_overlap")
            additive_schwartz_overlap = value;
        else if (property == "gmres_substeps")
            solver_gmres_substeps = value;
        else std::cout << "Parameter " << property << " of integral type is unknown" << std::endl;
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
