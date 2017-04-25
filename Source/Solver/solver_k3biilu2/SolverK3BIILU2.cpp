#include "SolverK3BIILU2.h"
#include "solver_k3biilu2.h"
#include <Source/Misc/utils.h>

namespace INMOST {

    SolverK3BIILU2::SolverK3BIILU2() {
        solver_data = NULL;
        matrix_data = NULL;
    }

    SolverInterface *SolverK3BIILU2::Copy(const SolverInterface *other) {
        if (other == NULL) {
            throw INMOST::SolverCopyNullException;
        }
        const SolverK3BIILU2 *k3other;
        try {
            k3other = dynamic_cast<const SolverK3BIILU2 *>(other);
        } catch (...) {
            throw INMOST::SolverCopyException;
        }
        SolverCopyDataK3biilu2(&solver_data, k3other->solver_data, k3other->communicator);
        if (k3other->matrix_data != NULL) {
            MatrixCopyDataK3biilu2(&matrix_data, k3other->matrix_data);
            SolverSetMatrixK3biilu2(solver_data, matrix_data, false, false);
        }
        return this;
    }

    void SolverK3BIILU2::Assign(const SolverInterface *other) {
        if (other == NULL) {
            throw INMOST::SolverAssignNullException;
        }
        const SolverK3BIILU2 *k3other;
        try {
            k3other = dynamic_cast<const SolverK3BIILU2 *>(other);
        } catch (...) {
            throw INMOST::SolverAssignException;
        }
        SolverAssignDataK3biilu2(solver_data, k3other->solver_data);
        if (k3other->matrix_data != NULL) {
            if (matrix_data != NULL) {
                MatrixAssignDataK3biilu2(matrix_data, k3other->matrix_data);
            } else {
                MatrixCopyDataK3biilu2(&matrix_data, k3other->matrix_data);
            }
            SolverSetMatrixK3biilu2(solver_data, matrix_data, false, false);
        }
    }

    void SolverK3BIILU2::Setup(int *argc, char ***argv, SolverParameters &p) {
        SolverInitDataK3biilu2(&solver_data, communicator, p.solverPrefix.c_str());
        if (p.internalFile != "") {
            SolverInitializeK3biilu2(solver_data, argc, argv, p.internalFile.c_str());
        } else {
            for (parameters_iterator_t parameter = p.parameters.begin(); parameter < p.parameters.end(); parameter++) {
                this->SetParameter((*parameter).first, (*parameter).second);
            }
        }
    }

    void SolverK3BIILU2::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) {
        bool modified_pattern = ModifiedPattern;
        //~ if( A.comm != comm ) throw DifferentCommunicatorInSolver;
        if (matrix_data == NULL) {
            MatrixInitDataK3biilu2(&matrix_data, A.GetCommunicator(), A.GetName().c_str());
            modified_pattern = true;
        }
        if (modified_pattern) {
            global_size = local_size = A.Size();

            MatrixDestroyDataK3biilu2(&matrix_data);
            MatrixInitDataK3biilu2(&matrix_data, A.GetCommunicator(), A.GetName().c_str());
            INMOST_DATA_ENUM_TYPE nnz = 0, k = 0, q = 1, shift = 0;
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
            for (Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++) {
                nnz += it->Size();
            }
            int *ja = (int *) malloc(sizeof(int) * nnz);
            double *values = (double *) malloc(sizeof(double) * nnz);
            ia[0] = shift;
            for (Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++) {
                for (Sparse::Row::iterator jt = it->Begin(); jt != it->End(); jt++) {
                    ja[k] = jt->first + 0;
                    values[k] = jt->second;
                    //std::cout<<"# q="<<q<<" k="<<k<<" ja="<<ja[k]<<" a="<<values[k]<<std::endl;//db!
                    k++;
                }
                shift += it->Size();
                ia[q++] = shift;
            }
            MatrixFillK3biilu2(matrix_data, n, nproc, ibl, ia, ja, values);
        } else {
            INMOST_DATA_ENUM_TYPE nnz = 0, k = 0;
            for (Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++) {
                nnz += it->Size();
            }
            double *values = (double *) malloc(sizeof(double) * nnz);
            k = 0;
            for (Sparse::Matrix::iterator it = A.Begin(); it != A.End(); it++) {
                for (Sparse::Row::iterator jt = it->Begin(); jt != it->End(); jt++) {
                    values[k++] = jt->second;
                }
            }
            MatrixFillValuesK3biilu2(matrix_data, values);
        }
        MatrixFinalizeK3biilu2(matrix_data);
        SolverSetMatrixK3biilu2(solver_data, matrix_data, modified_pattern, OldPreconditioner);
    }

    bool SolverK3BIILU2::Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL) {
        INMOST_DATA_ENUM_TYPE vbeg, vend;
        RHS.GetInterval(vbeg, vend);
        vector_k3biilu2 *rhs_data = NULL;
        vector_k3biilu2 *solution_data = NULL;

        VectorInitDataK3biilu2(&rhs_data, RHS.GetCommunicator(), RHS.GetName().c_str());
        VectorPreallocateK3biilu2(rhs_data, local_size);

        VectorInitDataK3biilu2(&solution_data, SOL.GetCommunicator(), SOL.GetName().c_str());
        VectorPreallocateK3biilu2(solution_data, local_size);

        VectorFillK3biilu2(rhs_data, &RHS[vbeg]);
        VectorFinalizeK3biilu2(rhs_data);

        VectorFillK3biilu2(solution_data, &SOL[vbeg]);
        VectorFinalizeK3biilu2(solution_data);

        bool result = SolverSolveK3biilu2(solver_data, rhs_data, solution_data);
        if (result) VectorLoadK3biilu2(solution_data, &SOL[vbeg]);
        iter_time = solver_data->dstat[9];
        time_prec = solver_data->dstat[7];
        return result;
    }

    bool SolverK3BIILU2::Clear() {
        if (matrix_data != NULL) {
            MatrixDestroyDataK3biilu2(&matrix_data);
        }
        SolverDestroyDataK3biilu2(&solver_data);
        return true;
    }

    bool SolverK3BIILU2::isMatrixSet() {
        return matrix_data != NULL;
    }

    std::string SolverK3BIILU2::GetParameter(std::string name) const {
        if (name == "time_prec") return INMOST::to_string(time_prec);
        if (name == "time_iter") return INMOST::to_string(iter_time);
        if (name == "time_total") return INMOST::to_string(time_prec + iter_time);
        if (name == "prec_density") return INMOST::to_string(solver_data->dstat[0]);
        if (name == "pivot_mod") return INMOST::to_string(solver_data->istat[0]);
#if !defined(SILENCE_SET_PARAMETER)
        std::cout << "Parameter " << name << " is unknown" << std::endl;
#endif
        return "";
    }

    void SolverK3BIILU2::SetParameter(std::string name, std::string value) {
        const char *val = value.c_str();
        if (name == "msglev") solver_data->pParIter->msglev = atoi(val);
        else if (name == "kovl") solver_data->pParams->ncycle = atoi(val);
        else if (name == "tau") {
            solver_data->pParams->tau1 = atof(val);
            solver_data->pParams->tau2 = -1.0;
        } else if (name == "nit") solver_data->pParIter->maxit = atoi(val);
        else if (name == "eps") solver_data->pParIter->eps = atof(val);
#if !defined(SILENCE_SET_PARAMETER)
        else std::cout << "Parameter " << name << " is unknown" << std::endl;
#endif
    }

    INMOST_DATA_ENUM_TYPE SolverK3BIILU2::Iterations() const {
        return static_cast<INMOST_DATA_ENUM_TYPE>(solver_data->istat[2]);
    }

    INMOST_DATA_REAL_TYPE SolverK3BIILU2::Residual() const {
        return solver_data->dstat[2];
    }

    const std::string SolverK3BIILU2::ReturnReason() const {
        return "Unspecified for K3BIILU2";
    }

    const std::string SolverK3BIILU2::SolverName() const {
        return "k3biilu2";
    }

    void SolverK3BIILU2::Finalize() {
        SolverFinalizeK3biilu2();
    }

    SolverK3BIILU2::~SolverK3BIILU2() {
        this->Clear();
    }

}
