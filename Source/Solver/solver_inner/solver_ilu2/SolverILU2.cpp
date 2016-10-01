#include "SolverILU2.h"

namespace INMOST {

    SolverILU2::SolverILU2() {
        Method *preconditioner = new ILU2_preconditioner(info);
        solver = new KSOLVER(preconditioner, info);
        matrix = NULL;
    }

    SolverILU2::SolverILU2(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverILU2::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) {
        if (matrix != NULL) {
            delete matrix;
        }
        matrix = new Sparse::Matrix(A);
        info.PrepareMatrix(*matrix, parameters.GetParameterEnum("additive_schwartz_overlap"));
        solver->ReplaceMAT(*matrix);

        solver->RealParameter(":tau") = parameters.GetParameterReal("drop_tolerance");
        solver->RealParameter(":tau2") = parameters.GetParameterReal("reuse_tolerance");
        solver->EnumParameter(":scale_iters") = parameters.GetParameterEnum("rescale_iterations");
        solver->EnumParameter(":fill") = static_cast<INMOST_DATA_ENUM_TYPE>(parameters.GetParameterReal("fill_level"));

        if (sizeof(KSOLVER) == sizeof(BCGSL_solver)) {
            solver->EnumParameter("levels") = parameters.GetParameterEnum("gmres_substeps");
        }

        if (!solver->isInitialized()) {
            solver->Initialize();
        }
    }

    const std::string SolverILU2::SolverName() const {
        return "inner_ilu2";
    }

    SolverILU2::~SolverILU2() {

    }

}