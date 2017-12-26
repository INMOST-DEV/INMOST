#include "SolverDDPQILUC2.h"

namespace INMOST {

    SolverDDPQILUC2::SolverDDPQILUC2() {
        Method *preconditioner = new ILUC_preconditioner(info);
        solver = new KSOLVER(preconditioner, info);
        matrix = NULL;
        rescale_iterations = 6;
        condition_estimation = 1;
        adapt_ddpq_tolerance = 1;
        schwartz_overlap = 1;
        gmres_substeps = 2;
        reorder_nnz = 1;
        drop_tolerance = 0.005;
        reuse_tolerance = 0.00005;
        ddpq_tolerance = 0.75;
    }

    SolverDDPQILUC2::SolverDDPQILUC2(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
        (void) other;
    }

    void SolverDDPQILUC2::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) {
        if (matrix != NULL) {
            delete matrix;
        }
        matrix = new Sparse::Matrix(A);
        info.PrepareMatrix(*matrix, schwartz_overlap);
        solver->ReplaceMAT(*matrix);

        solver->RealParameter(":tau") = drop_tolerance;
        solver->RealParameter(":tau2") = reuse_tolerance;
        solver->EnumParameter(":scale_iters") = rescale_iterations;
        solver->RealParameter(":ddpq_tau") = ddpq_tolerance;
        solver->EnumParameter(":reorder_nnz") = reorder_nnz;
        solver->EnumParameter(":estimator") = condition_estimation;
        solver->EnumParameter(":ddpq_tau_adapt") = adapt_ddpq_tolerance;

        if (sizeof(KSOLVER) == sizeof(BCGSL_solver)) {
            solver->EnumParameter("levels") = gmres_substeps;
        }

        if (!solver->isInitialized()) {
            solver->Initialize();
        }
    }

    void SolverDDPQILUC2::SetParameter(std::string name, std::string value) {
        const char *val = value.c_str();
        if (name == "rescale_iterations") rescale_iterations = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "condition_estimation") condition_estimation = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "adapt_ddpq_tolerance") adapt_ddpq_tolerance = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "schwartz_overlap") schwartz_overlap = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "gmres_substeps") gmres_substeps = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "reorder_nonzeros") reorder_nnz = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "drop_tolerance") drop_tolerance = atof(val);
        else if (name == "reuse_tolerance") reuse_tolerance = atof(val);
        else if (name == "ddpq_tolerance") ddpq_tolerance = atof(val);
        else SolverInner::SetParameter(name, value);
    }

    const std::string SolverDDPQILUC2::SolverName() const {
        return "inner_ddpqiluc2";
    }

    SolverDDPQILUC2::~SolverDDPQILUC2() {

    }

}
