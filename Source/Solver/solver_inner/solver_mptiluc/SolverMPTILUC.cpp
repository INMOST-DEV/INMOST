#include "SolverMPTILUC.h"

namespace INMOST {

    SolverMPTILUC::SolverMPTILUC() {
        Method *preconditioner = new MTILUC_preconditioner(info);
        solver = new KSOLVER(preconditioner, info);
        matrix = NULL;
    }

    SolverMPTILUC::SolverMPTILUC(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
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

    const std::string SolverMPTILUC::SolverName() const {
        return "inner_mptiluc";
    }

    SolverMPTILUC::~SolverMPTILUC() {

    }

}
