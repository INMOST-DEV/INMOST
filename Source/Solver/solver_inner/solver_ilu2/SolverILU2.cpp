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
        info.PrepareMatrix(*matrix, overlap);
        solver->ReplaceMAT(*matrix);

        solver->RealParameter(":tau") = tau;
        solver->RealParameter(":tau2") = tau2;
        solver->EnumParameter(":scale_iters") = scale_iters;
        solver->EnumParameter(":fill") = fill;

        if (sizeof(KSOLVER) == sizeof(BCGSL_solver)) {
            solver->EnumParameter("levels") = ell;
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