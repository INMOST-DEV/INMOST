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
        info.PrepareMatrix(*matrix, overlap);
        solver->ReplaceMAT(*matrix);

        solver->RealParameter(":tau") = tau;
        solver->RealParameter(":tau2") = tau2;
        solver->EnumParameter(":scale_iters") = scale_iters;
        solver->EnumParameter(":estimator") = condest;

        if (sizeof(KSOLVER) == sizeof(BCGSL_solver)) {
            solver->EnumParameter("levels") = ell;
        }

        if (!solver->isInitialized()) {
            solver->Initialize();
        }
    }

    bool SolverMPTILUC::Solve(Sparse::Vector &RHS, Sparse::Vector &SOL) {
        solver->EnumParameter("maxits") = iters;
        solver->RealParameter("rtol") = rtol;
        solver->RealParameter("atol") = atol;
        solver->RealParameter("divtol") = dtol;

        bool solve = solver->Solve(RHS, SOL);
        return solve;
    }

    const std::string SolverMPTILUC::SolverName() const {
        return "inner_mptiluc";
    }

    SolverMPTILUC::~SolverMPTILUC() {

    }

}
