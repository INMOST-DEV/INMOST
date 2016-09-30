#include "SolverAztec.h"

namespace INMOST {

    SolverAztec::SolverAztec() {
        problem = NULL;
        matrix = NULL;
    }

    SolverAztec::SolverAztec(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverAztec::Assign(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverAztec::Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix) {
        problem = new std::pair<std::string, Epetra_LinearProblem>;
        problem->first = prefix;
    }

    void SolverAztec::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) {
#if defined(USE_MPI)
        Epetra_MpiComm Comm(A.GetCommunicator());
#else
        Epetra_SerialComm Comm;
#endif
        bool modified_pattern = ModifiedPattern;
        INMOST_DATA_ENUM_TYPE mbeg, mend;
        A.GetInterval(mbeg, mend);
        bool refill = true;
        if (matrix == NULL || modified_pattern) {
            if (matrix != NULL) {
                delete matrix;
                matrix = NULL;
            }
            local_size = A.Size();
#if defined(USE_MPI)
            MPI_Allreduce(&local_size, &global_size, 1, INMOST_MPI_DATA_ENUM_TYPE, MPI_SUM, communicator);
#else
            global_size = local_size;
#endif
            int k;
            int imbeg = static_cast<int>(mbeg);
            int imend = static_cast<int>(mend);
            int *temporary = new int[imend - imbeg];
            for (k = imbeg; k < imend; ++k) {
                temporary[k - imbeg] = k;
            }
            Epetra_Map Map(static_cast<int>(global_size), static_cast<int>(local_size), temporary, 0, Comm);
            k = 0;
            for (k = imbeg; k < imend; ++k) {
                temporary[k - imbeg] = A[k].Size();
            }
            matrix = new Epetra_CrsMatrix(Copy, Map, temporary, true);
            delete[] temporary;
            refill = false;
        }
        {
            assert(matrix);
            unsigned max = 0;
            for (Sparse::Matrix::iterator it = A.Begin(); it != A.End(); ++it)
                if (it->Size() > max) max = it->Size();
            int *col_positions = new int[max];
            double *col_values = new double[max];
            for (INMOST_DATA_ENUM_TYPE k = mbeg; k < mend; ++k) {
                INMOST_DATA_ENUM_TYPE m = 0;
                for (INMOST_DATA_ENUM_TYPE i = 0; i < A[k].Size(); i++) {
                    col_positions[m] = A[k].GetIndex(i);
                    col_values[m] = A[k].GetValue(i);
                    m++;
                }
                int ret;
                if (refill) {
                    ret = matrix->ReplaceGlobalValues(k, m, col_values, col_positions);
                } else {
                    ret = matrix->InsertGlobalValues(k, m, col_values, col_positions);
                }
                if (ret != 0)
                    std::cout << "Problem reported by Trilinos! return code " << ret << std::endl;
            }
            delete[] col_positions;
            delete[] col_values;
            matrix->FillComplete();
            assert(problem);
            problem->second.SetOperator(matrix);
        }
    }

    bool SolverAztec::Solve(INMOST::Sparse::Vector &RHS, INMOST::Sparse::Vector &SOL) {

    }

    bool SolverAztec::Clear() {
        if (problem != NULL) {
            delete problem;
            problem = NULL;
        }
        if (matrix != NULL) {
            delete matrix;
            matrix = NULL;
        }
        return true;
    }

    bool SolverAztec::isMatrixSet() {
        return matrix != NULL;
    }

    INMOST_DATA_REAL_TYPE SolverAztec::GetPropertyReal(std::string property) const {

    }

    INMOST_DATA_ENUM_TYPE SolverAztec::GetPropertyEnum(std::string property) const {

    }

    void SolverAztec::SetPropertyReal(std::string property, INMOST_DATA_REAL_TYPE value) {

    }

    void SolverAztec::SetPropertyEnum(std::string property, INMOST_DATA_ENUM_TYPE value) {

    }

    const INMOST_DATA_ENUM_TYPE SolverAztec::Iterations() const {

    }

    const INMOST_DATA_REAL_TYPE SolverAztec::Residual() const {

    }

    const std::string SolverAztec::ReturnReason() const {

    }

    const std::string SolverAztec::SolverName() const {

    }

    void SolverAztec::Finalize() {

    }

    SolverAztec::~SolverAztec() {
        this->Clear();
    }
}