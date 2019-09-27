#include "SolverMLMPTILUC.h"

namespace INMOST
{

    SolverMLMPTILUC::SolverMLMPTILUC()
	{
        Method *preconditioner = new MLMTILUC_preconditioner(info);
        solver = new KSOLVER(preconditioner, info);
        matrix = NULL;
        rescale_iterations = 6;
        condition_estimation = 1;
        schwartz_overlap = 1;
        gmres_substeps = 2;
        reorder_nnz = 1;

        drop_tolerance = 0.005;
        reuse_tolerance = 0.00005;
		
		pivot_condition = 1.0e+6;
		pivot_diag = 1.0e+6;
        fill_level = 3;
        verbosity = 0;
    }

    SolverMLMPTILUC::SolverMLMPTILUC(const SolverInterface *other)
	{
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
        (void)other;
    }

    void SolverMLMPTILUC::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner)
	{
        if (matrix != NULL)
            delete matrix;
        matrix = new Sparse::Matrix(A);
        info.PrepareMatrix(*matrix, schwartz_overlap);
        solver->ReplaceMAT(*matrix);

        solver->RealParameter(":tau") = drop_tolerance;
        solver->RealParameter(":tau2") = reuse_tolerance;
        solver->EnumParameter(":scale_iters") = rescale_iterations;
        solver->EnumParameter(":estimator") = condition_estimation;
        solver->EnumParameter(":verbosity") = verbosity;
        solver->EnumParameter("verbosity") = verbosity;
		solver->RealParameter(":pivot_cond") = pivot_condition;
		solver->RealParameter(":pivot_diag") = pivot_diag;

        if (sizeof(KSOLVER) == sizeof(BCGSL_solver))
		    solver->EnumParameter("levels") = gmres_substeps;
		
        if (!solver->isInitialized())
            solver->Initialize();

        (void)ModifiedPattern;
        (void)OldPreconditioner;
    }

    void SolverMLMPTILUC::SetParameter(std::string name, std::string value)
	{
        const char *val = value.c_str();
        if (name == "rescale_iterations") rescale_iterations = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "condition_estimation") condition_estimation = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "schwartz_overlap") schwartz_overlap = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "gmres_substeps") gmres_substeps = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "reorder_nonzeros") reorder_nnz = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "fill_level") fill_level = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "drop_tolerance") drop_tolerance = atof(val);
        else if (name == "reuse_tolerance") reuse_tolerance = atof(val);
		else if (name == "pivot_condition") pivot_condition = atof(val);
		else if (name == "pivot_diag") pivot_diag = atof(val);
		else if (name == "verbosity") verbosity = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else SolverInner::SetParameter(name, value);
    }

    const std::string SolverMLMPTILUC::SolverName() const
	{
        return "inner_mptiluc";
    }

    SolverMLMPTILUC::~SolverMLMPTILUC() {}
}
