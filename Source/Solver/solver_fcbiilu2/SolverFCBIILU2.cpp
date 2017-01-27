#include <Source/Misc/utils.h>
#include "SolverFCBIILU2.h"
#include "solver_fcbiilu2.h"

namespace INMOST {

    SolverFCBIILU2::SolverFCBIILU2() {
        solver_data = NULL;
        matrix_data = NULL;
    }

    SolverInterface *SolverFCBIILU2::Copy(const SolverInterface *other) {
        if (other == NULL) {
            throw INMOST::SolverCopyNullException;
        }
        const SolverFCBIILU2 *fcother;
        try {
            fcother = dynamic_cast<const SolverFCBIILU2 *>(other);
        } catch (...) {
            throw INMOST::SolverCopyException;
        }
        SolverCopyDataFcbiilu2(&solver_data, fcother->solver_data, communicator);
        if (fcother->matrix_data != NULL) {
            MatrixCopyDataFcbiilu2(&matrix_data, fcother->matrix_data);
            SolverSetMatrixFcbiilu2(solver_data, matrix_data, false, false);
        }
        return this;
    }

    void SolverFCBIILU2::Assign(const SolverInterface *other) {
        if (other == NULL) {
            throw INMOST::SolverAssignNullException;
        }
        const SolverFCBIILU2 *fcother;
        try {
            fcother = dynamic_cast<const SolverFCBIILU2 *>(other);
        } catch (...) {
            throw INMOST::SolverAssignException;
        }
        SolverAssignDataFcbiilu2(solver_data, fcother->solver_data);
        if (fcother->matrix_data != NULL) {
            if (matrix_data != NULL) {
                MatrixAssignDataFcbiilu2(matrix_data, fcother->matrix_data);
            } else {
                MatrixCopyDataFcbiilu2(&matrix_data, fcother->matrix_data);
            }
            SolverSetMatrixFcbiilu2(solver_data, matrix_data, false, false);
        }
    }

    void SolverFCBIILU2::Setup(int *argc, char ***argv, SolverParameters &p) {
        SolverInitDataFcbiilu2(&solver_data, communicator, p.solverPrefix.c_str());
        solver_data->kovl = 0;    // number of overlap layers: kovl=0,1,2,...
        solver_data->tau = 3e-3;  // the ILU2 precision (for the submatrix factorization); tau=3e-3
        solver_data->eps = 1e-5;  // the residual precision: ||r|| < eps * ||b||; eps=1e-6
        solver_data->nit = 999;   // number of iterations permitted; nit=999
        solver_data->msglev = 2;  // messages level; msglev=0 for silent; msglev=1 to output solution statistics
        SolverInitializeFcbiilu2(solver_data, argc, argv, p.internalFile.c_str());
        for (parameters_iterator_t parameter = p.parameters.begin(); parameter < p.parameters.end(); parameter++) {
            this->SetParameter((*parameter).first, (*parameter).second);
        }
    }

    void SolverFCBIILU2::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) {
        bool modified_pattern = ModifiedPattern;
        //~ if( A.comm != comm ) throw DifferentCommunicatorInSolver;
        if (matrix_data == NULL) {
            MatrixInitDataFcbiilu2(&matrix_data, A.GetCommunicator(), A.GetName().c_str());
            modified_pattern = true;
        }
        if (modified_pattern) {
            global_size = local_size = A.Size();

            MatrixDestroyDataFcbiilu2(&matrix_data);
            MatrixInitDataFcbiilu2(&matrix_data, A.GetCommunicator(), A.GetName().c_str());
            INMOST_DATA_ENUM_TYPE nnz = 0, k = 0, q = 1, shift = 1;
            int nproc, myid;
#if defined(USE_MPI)
            MPI_Comm_size(A.GetCommunicator(), &nproc);
            MPI_Comm_rank(A.GetCommunicator(), &myid);
#else
            nproc = 1;
            myid = 0;
#endif
            int *ibl = (int *) malloc(sizeof(int) * (nproc + 1));
            ibl[0] = 0;
            int n = A.Size();
#if defined(USE_MPI)
            INMOST_DATA_ENUM_TYPE mbeg, mend;
            A.GetInterval(mbeg, mend);
            n = mend - mbeg;
            int block_end = mend;
            MPI_Allgather(&block_end, 1, MPI_INT, &ibl[1], 1, MPI_INT, A.GetCommunicator());
#else
            ibl[1] = n;
#endif
            int *ia = (int *) malloc(sizeof(int) * (n + 1));
            for (Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++) nnz += it->Size();
            int *ja = (int *) malloc(sizeof(int) * nnz);
            double *values = (double *) malloc(sizeof(double) * nnz);
            ia[0] = shift;
            for (Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++) {
                for (Sparse::Row::iterator jt = it->Begin(); jt != it->End(); jt++) {
                    ja[k] = jt->first + 1;
                    values[k] = jt->second;
                    //std::cout<<"# q="<<q<<" k="<<k<<" ja="<<ja[k]<<" a="<<values[k]<<std::endl;//db!
                    k++;
                }
                shift += it->Size();
                ia[q++] = shift;
            }
            MatrixFillFcbiilu2(matrix_data, n, nproc, ibl, ia, ja, values);
        } else {
            INMOST_DATA_ENUM_TYPE nnz = 0, k = 0;
            for (Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++) nnz += it->Size();
            double *values = (double *) malloc(sizeof(double) * nnz);
            k = 0;
            for (Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++)
                for (Sparse::Row::iterator jt = it->Begin(); jt != it->End(); jt++)
                    values[k++] = jt->second;
            MatrixFillValuesFcbiilu2(matrix_data, values);
        }
        MatrixFinalizeFcbiilu2(matrix_data);
        SolverSetMatrixFcbiilu2(solver_data, matrix_data, modified_pattern, OldPreconditioner);
        time_prec = solver_data->dstat[7];
    }

    bool SolverFCBIILU2::Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL) {
        INMOST_DATA_ENUM_TYPE vbeg, vend;
        RHS.GetInterval(vbeg, vend);

        vector_fcbiilu2 *rhs_data = NULL;
        VectorInitDataFcbiilu2(&rhs_data, RHS.GetCommunicator(), RHS.GetName().c_str());
        VectorPreallocateFcbiilu2(rhs_data, local_size);

        vector_fcbiilu2 *solution_data = NULL;
        VectorInitDataFcbiilu2(&solution_data, SOL.GetCommunicator(), SOL.GetName().c_str());
        VectorPreallocateFcbiilu2(solution_data, local_size);
        VectorFillFcbiilu2(rhs_data, &RHS[vbeg]);
        VectorFinalizeFcbiilu2(rhs_data);

        VectorFillFcbiilu2(solution_data, &SOL[vbeg]);
        VectorFinalizeFcbiilu2(solution_data);
        bool result = SolverSolveFcbiilu2(solver_data, rhs_data, solution_data);
        if (result) VectorLoadFcbiilu2(solution_data, &SOL[vbeg]);
        iter_time = solver_data->dstat[9];
        return result;
    }

    bool SolverFCBIILU2::Clear() {
        if (matrix_data != NULL) {
            MatrixDestroyDataFcbiilu2(&matrix_data);
        }
        SolverDestroyDataFcbiilu2(&solver_data);
        return true;
    }

    bool SolverFCBIILU2::isMatrixSet() {
        return matrix_data != NULL;
    }

    std::string SolverFCBIILU2::GetParameter(std::string name) const {
        if (name == "time_prec") return INMOST::to_string(time_prec);
        if (name == "time_iter") return INMOST::to_string(iter_time);
        if (name == "time_total") return INMOST::to_string(time_prec + iter_time);
        if (name == "prec_density") return INMOST::to_string(solver_data->dstat[0]);
        if (name == "pivot_mod") return INMOST::to_string(solver_data->istat[0]);
        std::cout << "Parameter " << name << " is unknown" << std::endl;
        return "";
    }

    void SolverFCBIILU2::SetParameter(std::string name, std::string value) {
        const char *val = value.c_str();
        if (name == "kovl") solver_data->kovl = atoi(val);
        else if (name == "tau") solver_data->tau = atof(val);
        else if (name == "eps") solver_data->eps = atof(val);
        else if (name == "nit") solver_data->nit = atoi(val);
        else if (name == "msglev") solver_data->msglev = atoi(val);
        else std::cout << "Parameter " << name << " is unknown" << std::endl;
    }

    const INMOST_DATA_ENUM_TYPE SolverFCBIILU2::Iterations() const {
        return static_cast<INMOST_DATA_ENUM_TYPE>(solver_data->ITER);
    }

    const INMOST_DATA_REAL_TYPE SolverFCBIILU2::Residual() const {
        return solver_data->RESID;
    }

    const std::string SolverFCBIILU2::ReturnReason() const {
        return "Unspecified for FCBIILU2";
    }

    const std::string SolverFCBIILU2::SolverName() const {
        return "fcbiilu2";
    }

    void SolverFCBIILU2::Finalize() {
        SolverFinalizeFcbiilu2();
    }

    SolverFCBIILU2::~SolverFCBIILU2() {
        this->Clear();
    }

}