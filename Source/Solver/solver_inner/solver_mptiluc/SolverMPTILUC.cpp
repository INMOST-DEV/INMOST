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
        info.PrepareMatrix(*matrix, parameters.GetParameter("additive_schwartz_overlap").unsigned_integer());
        solver->ReplaceMAT(*matrix);

        solver->RealParameter(":tau") = parameters.GetParameter("drop_tolerance").real();
        solver->RealParameter(":tau2") = parameters.GetParameter("reuse_tolerance").real();
        solver->EnumParameter(":scale_iters") = parameters.GetParameter("rescale_iterations").unsigned_integer();
        solver->EnumParameter(":estimator") = parameters.GetParameter("condition_estimation").unsigned_integer();

        if (sizeof(KSOLVER) == sizeof(BCGSL_solver)) {
            solver->EnumParameter("levels") = parameters.GetParameter("gmres_substeps").unsigned_integer();
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
