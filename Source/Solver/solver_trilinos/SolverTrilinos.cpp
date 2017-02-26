#include "SolverTrilinos.h"

namespace INMOST {

    SolverTrilinos::SolverTrilinos() {
        maximum_iterations = 2500;
        schwartz_overlap = 1;

        rtol = 1.0e-12;
        drop_tolerance = 0.005;
        fill_level = 3;

        Epetra_problem = NULL;
        matrix = NULL;
    }

    SolverInterface *SolverTrilinos::Copy(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverTrilinos::Assign(const SolverInterface *other) {
        //You should not really want to copy solver's information
        throw INMOST::SolverUnsupportedOperation;
    }

    void SolverTrilinos::Setup(int *argc, char ***argv, SolverParameters &p) {
        this->Epetra_problem = new std::pair<std::string, Epetra_LinearProblem>;
        this->Epetra_problem->first = p.solverPrefix;
        this->parameters_file = p.internalFile;
        for (parameters_iterator_t parameter = p.parameters.begin(); parameter < p.parameters.end(); parameter++) {
            this->SetParameter((*parameter).first, (*parameter).second);
        }
    }

    void SolverTrilinos::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) {
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
            matrix = new Epetra_CrsMatrix(::Copy, Map, temporary, true);
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
            assert(Epetra_problem);
            Epetra_problem->second.SetOperator(matrix);
        }
    }

    bool SolverTrilinos::Clear() {
        if (Epetra_problem != NULL) {
            delete Epetra_problem;
            Epetra_problem = NULL;
        }
        if (matrix != NULL) {
            delete matrix;
            matrix = NULL;
        }
        return true;
    }

    bool SolverTrilinos::isMatrixSet() {
        return matrix != NULL;
    }

    std::string SolverTrilinos::GetParameter(std::string name) const {
        if(name == "maximum_iterations" ) return to_string(maximum_iterations);
        else if( name == "schwartz_overlap" ) return to_string(schwartz_overlap);
        else if( name == "relative_tolerance") return to_string(rtol);
        else if( name == "drop_tolerance") return to_string(drop_tolerance);
        else if( name == "fill_level") return to_string(fill_level);
        else {
#if !defined(SILENCE_SET_PARAMETER)
            std::cout << "Parameter " << name << " is unknown (Use internal file for all parameters)" << std::endl;
#endif
            return "";
        }
    }

    void SolverTrilinos::SetParameter(std::string name, std::string value) {
        const char *val = value.c_str();
        if(name == "maximum_iterations" ) maximum_iterations = atoi(val);
        else if( name == "schwartz_overlap" ) schwartz_overlap = atoi(val);
        else if( name == "relative_tolerance") rtol = atof(val);
        else if( name == "drop_tolerance") drop_tolerance = atof(val);
        else if( name == "fill_level") fill_level = atoi(val);
#if !defined(SILENCE_SET_PARAMETER)
        else std::cout << "Parameter " << name << " is unknown (Use internal file for all parameters)" << std::endl;
#endif
    }

    const INMOST_DATA_ENUM_TYPE SolverTrilinos::Iterations() const {
        return lastIterations;
    }

    const INMOST_DATA_REAL_TYPE SolverTrilinos::Residual() const {
        return lastResidual;
    }

    const std::string SolverTrilinos::ReturnReason() const {
        return returnReason;
    }

    void SolverTrilinos::Finalize() {

    }

    SolverTrilinos::~SolverTrilinos() {
        this->Clear();
    }

    void SolverTrilinos::TrilinosCheckStatus(int status_id, bool &success, std::string &reason) {
        success = true;
        switch (status_id) {
            case AZ_normal:
                reason = "User requested convergence criteria is satisfied.";
                break;
            case AZ_param:
                reason = "User requested option is not availible.";
                success = false;
                break;
            case AZ_breakdown:
                reason = "Numerical breakdown occurred.";
                success = false;
                break;
            case AZ_loss:
                reason = "Numerical loss precision occurred.";
                success = false;
                break;
            case AZ_ill_cond:
                reason = "The Hessenberg matrix within GMRES is illconditioned. "
                        "This could be caused by a number "
                        "of reasons. For example, the preconditioning "
                        "matrix could be nearly singular due to an unstable "
                        "factorization (note: pivoting is not implemented "
                        "in any of the incomplete factorizations). "
                        "Ill-conditioned Hessenberg matrices could also "
                        "arise from a singular application matrix. In this "
                        "case, GMRES tries to compute a least-squares "
                        "solution.";
                success = false;
                break;
            case AZ_maxits:
                reason = "Maximum iterations taken without convergence.";
                success = false;
                break;
            default: {
                std::stringstream str;
                str << "reason code " << status_id
                    << " was not specified in manual by the time, reffer to Trilinos manual.";
                reason = str.str();
            }
                break;
        }
    }
}