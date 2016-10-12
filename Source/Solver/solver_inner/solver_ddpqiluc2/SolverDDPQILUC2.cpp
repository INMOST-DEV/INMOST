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
        info.PrepareMatrix(*matrix, parameters.get<INMOST_DATA_ENUM_TYPE>("additive_schwartz_overlap"));
        solver->ReplaceMAT(*matrix);

        solver->RealParameter(":tau") = parameters.get<INMOST_DATA_REAL_TYPE>("drop_tolerance");
        solver->RealParameter(":tau2") = parameters.get<INMOST_DATA_REAL_TYPE>("reuse_tolerance");
        solver->EnumParameter(":scale_iters") = parameters.get<INMOST_DATA_ENUM_TYPE>("rescale_iterations");
        solver->RealParameter(":ddpq_tau") = parameters.get<INMOST_DATA_REAL_TYPE>("ddpq_tolerance");
        solver->EnumParameter(":reorder_nnz") = parameters.get<INMOST_DATA_ENUM_TYPE>("reorder_nonzero");
        solver->EnumParameter(":estimator") = parameters.get<INMOST_DATA_ENUM_TYPE>("condition_estimation");
        solver->EnumParameter(":ddpq_tau_adapt") = parameters.get<INMOST_DATA_ENUM_TYPE>("adapt_ddpq_tolerance");

        if (sizeof(KSOLVER) == sizeof(BCGSL_solver)) {
            solver->EnumParameter("levels") = parameters.get<INMOST_DATA_ENUM_TYPE>("gmres_substeps");
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
