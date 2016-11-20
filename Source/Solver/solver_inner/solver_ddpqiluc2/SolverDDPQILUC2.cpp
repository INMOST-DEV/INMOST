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
        info.PrepareMatrix(*matrix, overlap);
        solver->ReplaceMAT(*matrix);

        solver->RealParameter(":tau") = tau;
        solver->RealParameter(":tau2") = tau2;
        solver->EnumParameter(":scale_iters") = scale_iters;
        solver->RealParameter(":ddpq_tau") = ddpqtol;
        solver->EnumParameter(":reorder_nnz") = reorder_nnz;
        solver->EnumParameter(":estimator") = condest;
        solver->EnumParameter(":ddpq_tau_adapt") = ddpqatol;

        if (sizeof(KSOLVER) == sizeof(BCGSL_solver)) {
            solver->EnumParameter("levels") = ell;
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
