#include "SolverILU2.h"

namespace INMOST {

    SolverILU2::SolverILU2() {
        additive_schwartz_overlap = DEFAULT_ADDITIVE_SCHWARTZ_OVERLAP;
        maximum_iterations = DEFAULT_MAXIMUM_ITERATIONS;
        absolute_tolerance = DEFAULT_ABSOLUTE_TOLERANCE;
        relative_tolerance = DEFAULT_RELATIVE_TOLERANCE;
        divergence_tolerance = DEFAULT_DIVERGENCE_TOLERANCE;

        preconditioner_drop_tolerance = DEFAULT_PRECONDITIONER_DROP_TOLERANCE;
        preconditioner_reuse_tolerance = DEFAULT_PRECONDITIONER_REUSE_TOLERANCE;
        preconditioner_rescale_iterations = DEFAULT_PRECONDITIONER_RESCALE_ITERS;
        preconditioner_fill_level = DEFAULT_PRECONDITIONER_FILL_LEVEL;

        Method *prec = new ILU2_preconditioner(info);
        solver = new BCGS_solver(prec, info);
        matrix = NULL;
    }

    SolverILU2::SolverILU2(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverILU2::Assign(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverILU2::Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix) {

    }

    void SolverILU2::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) {
        Sparse::Matrix *m = new Sparse::Matrix(A);
        info.PrepareMatrix(*m, additive_schwartz_overlap);
        solver->ReplaceMAT(*m);
        if (matrix != NULL) {
            delete matrix;
        }
        matrix = m;

        solver->RealParameter(":tau") = preconditioner_drop_tolerance;
        solver->RealParameter(":tau2") = preconditioner_reuse_tolerance;
        solver->EnumParameter(":scale_iters") = preconditioner_rescale_iterations;
        solver->EnumParameter(":fill") = static_cast<INMOST_DATA_ENUM_TYPE>(preconditioner_fill_level);

        if (!solver->isInitialized()) {
            solver->Initialize();
        }
    }

    bool SolverILU2::Solve(Sparse::Vector &RHS, Sparse::Vector &SOL) {
        solver->EnumParameter("maxits") = maximum_iterations;
        solver->RealParameter("rtol") = relative_tolerance;
        solver->RealParameter("atol") = absolute_tolerance;
        solver->RealParameter("divtol") = divergence_tolerance;

        return solver->Solve(RHS, SOL);
    }

    bool SolverILU2::Clear() {
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

    bool SolverILU2::isMatrixSet() {
        return matrix != NULL;
    }

    INMOST_DATA_REAL_TYPE SolverILU2::GetPropertyReal(std::string property) const {
        return solver->RealParameter(property);
    }

    INMOST_DATA_ENUM_TYPE SolverILU2::GetPropertyEnum(std::string property) const {
        return solver->EnumParameter(property);
    }

    void SolverILU2::SetPropertyReal(std::string property, INMOST_DATA_REAL_TYPE value) {
        //Maybe we should use explicit names?, like maxits ..etc
        //like solver->RealParameter(property) = value;
        if (property == "absolute_tolerance") {
            absolute_tolerance = value;
        } else if (property == "relative_tolerance") {
            relative_tolerance = value;
        } else if (property == "divergence_tolerance") {
            divergence_tolerance = value;
        } else if (property == "drop_tolerance") {
            preconditioner_drop_tolerance = value;
        } else if (property == "reuse_tolerance") {
            preconditioner_reuse_tolerance = value;
        } else if (property == "fill_level") {
            preconditioner_fill_level = value;
        } else {
            std::cout << "Parameter " << property << " of real type is unknown" << std::endl;
        }
    }

    void SolverILU2::SetPropertyEnum(std::string property, INMOST_DATA_ENUM_TYPE value) {
        //Maybe we should use explicit names?, like maxits ..etc
        //like solver->EnumParameter(property) = value;
        if (property == "maximum_iterations") {
            maximum_iterations = value;
        } else if (property == "rescale_iterations") {
            preconditioner_rescale_iterations = value;
        } else if (property == "schwartz_overlap") {
            additive_schwartz_overlap = value;
        } else std::cout << "Parameter " << property << " of integral type is unknown" << std::endl;
    }

    const INMOST_DATA_ENUM_TYPE SolverILU2::Iterations() const {
        return solver->GetIterations();
    }

    const INMOST_DATA_REAL_TYPE SolverILU2::Residual() const {
        return solver->GetResidual();
    }

    const std::string SolverILU2::ReturnReason() const {
        return solver->GetReason();
    }

    const std::string SolverILU2::SolverName() const {
        return "inner_ilu2";
    }

    void SolverILU2::Finalize() {

    }

    SolverILU2::~SolverILU2() {
        Clear();
    }

}