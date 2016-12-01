#include "SolverFCBIILU2.h"

namespace INMOST {

    SolverFCBIILU2::SolverFCBIILU2() {

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
        SolverInitializeFcbiilu2(argc, argv, p.internalFile.c_str());
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
    }

    bool SolverFCBIILU2::Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL) {
        INMOST_DATA_ENUM_TYPE vbeg, vend;
        RHS.GetInterval(vbeg, vend);

        void *rhs_data = NULL;
        VectorInitDataFcbiilu2(&rhs_data, RHS.GetCommunicator(), RHS.GetName().c_str());
        VectorPreallocateFcbiilu2(rhs_data, local_size);

        void *solution_data = NULL;
        VectorInitDataFcbiilu2(&solution_data, SOL.GetCommunicator(), SOL.GetName().c_str());
        VectorPreallocateFcbiilu2(solution_data, local_size);
        VectorFillFcbiilu2(rhs_data, &RHS[vbeg]);
        VectorFinalizeFcbiilu2(rhs_data);

        VectorFillFcbiilu2(solution_data, &SOL[vbeg]);
        VectorFinalizeFcbiilu2(solution_data);
        bool result = SolverSolveFcbiilu2(solver_data, rhs_data, solution_data);
        if (result) VectorLoadFcbiilu2(solution_data, &SOL[vbeg]);
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
        std::cout << "SolverFCBIILU2::GetParameter unsupported operation" << std::endl;
        //throw INMOST::SolverUnsupportedOperation;
        return "";
    }

    void SolverFCBIILU2::SetParameter(std::string name, std::string value) {
        std::cout << "SolverFCBIILU2::SetParameter unsupported operation" << std::endl;
        //throw INMOST::SolverUnsupportedOperation;
    }

    const INMOST_DATA_ENUM_TYPE SolverFCBIILU2::Iterations() const {
        return static_cast<INMOST_DATA_ENUM_TYPE>(SolverIterationNumberFcbiilu2(solver_data));
    }

    const INMOST_DATA_REAL_TYPE SolverFCBIILU2::Residual() const {
        return SolverResidualNormFcbiilu2(solver_data);
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