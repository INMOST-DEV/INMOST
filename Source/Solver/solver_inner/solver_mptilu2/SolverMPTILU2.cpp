#include "SolverMPTILU2.h"

namespace INMOST {

    SolverMPTILU2::SolverMPTILU2() {
        Method *preconditioner = new MTILU2_preconditioner(info);
        solver = new KSOLVER(preconditioner, info);
        matrix = NULL;
    }

    SolverMPTILU2::SolverMPTILU2(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverMPTILU2::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) {
        if (matrix != NULL) {
            delete matrix;
        }
        matrix = new Sparse::Matrix(A);
        info.PrepareMatrix(*matrix, overlap);
        solver->ReplaceMAT(*matrix);

        solver->RealParameter(":tau") = tau;
        solver->RealParameter(":tau2") = tau2;
        solver->EnumParameter(":scale_iters") = scale_iters;

        if (sizeof(KSOLVER) == sizeof(BCGSL_solver)) {
            solver->EnumParameter("levels") = ell;
        }

        if (!solver->isInitialized()) {
            solver->Initialize();
        }
    }

    const std::string SolverMPTILU2::SolverName() const {
        return "inner_mptilu2";
    }

    SolverMPTILU2::~SolverMPTILU2() {

    }

}
