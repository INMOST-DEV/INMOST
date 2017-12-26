#include "SolverILU2.h"

namespace INMOST {

    SolverILU2::SolverILU2() {
        Method *preconditioner = new ILU2_preconditioner(info);
        solver = new KSOLVER(preconditioner, info);
        matrix = NULL;
        rescale_iterations = 6;
        schwartz_overlap = 1;
        gmres_substeps = 2;
        reorder_nnz = 1;
        drop_tolerance = 0.005;
        reuse_tolerance = 0.00005;
        fill_level = 3;
    }

    SolverILU2::SolverILU2(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
        (void) other;
    }

    void SolverILU2::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) {
        if (matrix != NULL) {
            delete matrix;
        }
        matrix = new Sparse::Matrix(A);
        info.PrepareMatrix(*matrix, schwartz_overlap);
        solver->ReplaceMAT(*matrix);

        solver->RealParameter(":tau") = drop_tolerance;
        solver->RealParameter(":tau2") = reuse_tolerance;
        solver->EnumParameter(":scale_iters") = rescale_iterations;
        solver->EnumParameter(":fill") = fill_level;

        if (sizeof(KSOLVER) == sizeof(BCGSL_solver)) {
            solver->EnumParameter("levels") = gmres_substeps;
        }

        if (!solver->isInitialized()) {
            solver->Initialize();
        }
    }

    void SolverILU2::SetParameter(std::string name, std::string value) {
        const char *val = value.c_str();
        if (name == "rescale_iterations") rescale_iterations = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "schwartz_overlap") schwartz_overlap = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "gmres_substeps") gmres_substeps = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "reorder_nonzeros") reorder_nnz = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "fill_level") fill_level = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "drop_tolerance") drop_tolerance = atof(val);
        else if (name == "reuse_tolerance") reuse_tolerance = atof(val);
        else SolverInner::SetParameter(name, value);
    }

    const std::string SolverILU2::SolverName() const {
        return "inner_ilu2";
    }

    SolverILU2::~SolverILU2() {

    }

}
