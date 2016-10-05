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

    void SolverDDPQILUC2::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) {
        if (matrix != NULL) {
            delete matrix;
        }
        matrix = new Sparse::Matrix(A);
        info.PrepareMatrix(*matrix, parameters.GetParameter("additive_schwartz_overlap").unsigned_integer());
        solver->ReplaceMAT(*matrix);

        solver->RealParameter(":tau") = parameters.GetParameter("drop_tolerance").real();
        solver->RealParameter(":tau2") = parameters.GetParameter("reuse_tolerance").real();
        solver->EnumParameter(":scale_iters") = parameters.GetParameter("rescale_iterations").unsigned_integer();
        solver->RealParameter(":ddpq_tau") = parameters.GetParameter("ddpq_tolerance").real();
        solver->EnumParameter(":reorder_nnz") = parameters.GetParameter("reorder_nonzero").unsigned_integer();
        solver->EnumParameter(":estimator") = parameters.GetParameter("condition_estimation").unsigned_integer();
        solver->EnumParameter(":ddpq_tau_adapt") = parameters.GetParameter("adapt_ddpq_tolerance").unsigned_integer();

        if (sizeof(KSOLVER) == sizeof(BCGSL_solver)) {
            solver->EnumParameter("levels") = parameters.GetParameter("gmres_substeps").unsigned_integer();
        }

        if (!solver->isInitialized()) {
            solver->Initialize();
        }
    }

    const std::string SolverDDPQILUC2::SolverName() const {
        return "inner_ddpqiluc2";
    }

    SolverDDPQILUC2::~SolverDDPQILUC2() {

    }

}
