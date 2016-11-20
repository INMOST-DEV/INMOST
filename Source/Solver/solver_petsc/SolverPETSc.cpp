#include "SolverPETSc.h"

namespace INMOST {

    unsigned int SolverPETSc::petscSolversCount = 0;

    SolverPETSc::SolverPETSc() {
        petscSolversCount++;
        this->ksp = NULL;
        this->matrix = NULL;
    }

    SolverPETSc::SolverPETSc(const SolverInterface *other) {
        petscSolversCount++;
        const SolverPETSc *solver = static_cast<const SolverPETSc *>(other);
        this->ksp = NULL;
        this->matrix = NULL;

        SolverCopyDataPetsc(&ksp, solver->ksp, solver->communicator);
        if (solver->matrix != NULL) {
            MatrixCopyDataPetsc(&matrix, solver->matrix);
            SolverSetMatrixPetsc(ksp, matrix, false, false);
        }
    }

    void SolverPETSc::Assign(const SolverInterface *other) {
        const SolverPETSc *other_solver = static_cast<const SolverPETSc *>(other);
        this->parametersFile = other_solver->parametersFile;
        SolverAssignDataPetsc(ksp, other_solver->ksp);
        if (other_solver->matrix != NULL) {
            if (matrix != NULL) {
                MatrixAssignDataPetsc(matrix, other_solver->matrix);
            } else {
                MatrixCopyDataPetsc(&matrix, other_solver->matrix);
            }
            SolverSetMatrixPetsc(ksp, matrix, false, false);
        }
    }

    void SolverPETSc::Initialize(int *argc, char ***argv, const char *parameters_file, std::string prefix) {
        this->parametersFile = std::string(parameters_file);
        SolverInitializePetsc(argc, argv, parameters_file);
        SolverInitDataPetsc(&ksp, this->communicator, prefix.c_str());
    }

    void SolverPETSc::SetMatrix(Sparse::Matrix &A, bool ModifiedPattern, bool OldPreconditioner) {
        bool modified_pattern = ModifiedPattern;
        if (!isMatrixSet()) {
            MatrixInitDataPetsc(&matrix, A.GetCommunicator(), A.GetName().c_str());
            modified_pattern = true;
        }
        INMOST_DATA_ENUM_TYPE mbeg, mend;
        A.GetInterval(mbeg, mend);
        if (modified_pattern) {
            local_size = A.Size();
#if defined(USE_MPI)
            MPI_Allreduce(&local_size, &global_size, 1, INMOST_MPI_DATA_ENUM_TYPE, MPI_SUM, this->communicator);
#else
            global_size = local_size;
#endif
            int max = 0;
            {
                int *diag_nonzeroes = new int[local_size];
                int *off_diag_nonzeroes = new int[local_size];
                unsigned k = 0;

                for (Sparse::Matrix::iterator it = A.Begin(); it != A.End(); ++it) {
                    diag_nonzeroes[k] = off_diag_nonzeroes[k] = 0;
                    for (INMOST_DATA_ENUM_TYPE i = 0; i < it->Size(); i++) {
                        INMOST_DATA_ENUM_TYPE index = it->GetIndex(i);
                        if (index < mbeg || index > mend - 1) {
                            off_diag_nonzeroes[k]++;
                        } else {
                            diag_nonzeroes[k]++;
                        }
                    }
                    if (diag_nonzeroes[k] + off_diag_nonzeroes[k] > max) {
                        max = diag_nonzeroes[k] + off_diag_nonzeroes[k];
                    }
                    k++;
                }
                MatrixPreallocatePetsc(matrix, local_size, global_size, diag_nonzeroes, off_diag_nonzeroes);
                delete[] diag_nonzeroes;
                delete[] off_diag_nonzeroes;
            }
            if (max > 0) {
                int *col_positions = new int[max];
                double *col_values = new double[max];
                unsigned k = 0, m;
                for (Sparse::Matrix::iterator it = A.Begin(); it != A.End(); ++it) {
                    m = 0;
                    for (INMOST_DATA_ENUM_TYPE i = 0; i < it->Size(); i++) {
                        col_positions[m] = it->GetIndex(i);
                        col_values[m] = it->GetValue(i);
                        m++;
                    }
                    MatrixFillPetsc(matrix, mbeg + k, m, col_positions, col_values);
                    k++;
                }
                delete[] col_positions;
                delete[] col_values;
            }
        } else {
            unsigned max = 0;
            for (Sparse::Matrix::iterator it = A.Begin(); it != A.End(); ++it) {
                if (it->Size() > max) {
                    max = it->Size();
                }
            }
            int *col_positions = new int[max];
            double *col_values = new double[max];
            unsigned k = 0, m;
            for (Sparse::Matrix::iterator it = A.Begin(); it != A.End(); ++it) {
                m = 0;
                for (INMOST_DATA_ENUM_TYPE i = 0; i < it->Size(); i++) {
                    col_positions[m] = it->GetIndex(i);
                    col_values[m] = it->GetValue(i);
                    m++;
                }
                MatrixFillPetsc(matrix, mbeg + k, m, col_positions, col_values);
                k++;
            }
            delete[] col_positions;
            delete[] col_values;
        }
        MatrixFinalizePetsc(matrix);

        if (parametersFile == "") {
            SolverSetDropTolerancePetsc(ksp, DEFAULT_PRECONDITIONER_DROP_TOLERANCE);
            SolverSetFillLevelPetsc(ksp, DEFAULT_PRECONDITIONER_FILL_LEVEL);
            SolverSetOverlapPetsc(ksp, DEFAULT_ADDITIVE_SCHWARTZ_OVERLAP);
        }

        SolverSetMatrixPetsc(ksp, matrix, modified_pattern, OldPreconditioner);
    }

    bool SolverPETSc::Solve(Sparse::Vector &RHS, Sparse::Vector &SOL) {
        Vec *rhs = NULL;
        VectorInitDataPetsc(&rhs, RHS.GetCommunicator(), RHS.GetName().c_str());
        VectorPreallocatePetsc(rhs, local_size, global_size);
        Vec *solution = NULL;
        VectorInitDataPetsc(&solution, SOL.GetCommunicator(), SOL.GetName().c_str());
        VectorPreallocatePetsc(solution, local_size, global_size);
        INMOST_DATA_ENUM_TYPE vbeg, vend;
        RHS.GetInterval(vbeg, vend);
        int *positions = new int[local_size];
        double *values = new double[local_size];
        unsigned int k = 0;
        for (Sparse::Vector::iterator it = RHS.Begin(); it != RHS.End(); ++it) {
            positions[k] = vbeg + k;
            values[k] = *it;
            k++;
        }
        VectorFillPetsc(rhs, local_size, positions, values);
        VectorFinalizePetsc(rhs);

        k = 0;
        for (Sparse::Vector::iterator it = SOL.Begin(); it != SOL.End(); ++it) {
            values[k] = *it;
            k++;
        }
        VectorFillPetsc(solution, local_size, positions, values);
        VectorFinalizePetsc(solution);

        if (parametersFile == "") {
            SolverSetTolerancesPetsc(ksp,
                                     DEFAULT_RELATIVE_TOLERANCE,
                                     DEFAULT_ABSOLUTE_TOLERANCE,
                                     DEFAULT_DIVERGENCE_TOLERANCE,
                                     DEFAULT_MAXIMUM_ITERATIONS);
        }

        bool result = SolverSolvePetsc(ksp, rhs, solution);
        if (result) {
            VectorLoadPetsc(solution, local_size, positions, values);
            k = 0;
            for (Sparse::Vector::iterator it = SOL.Begin(); it != SOL.End(); ++it) {
                *it = values[k];
                k++;
            }
        }
        VectorDestroyDataPetsc(&rhs);
        VectorDestroyDataPetsc(&solution);
        delete[] positions;
        delete[] values;
        return result;
    }

    bool SolverPETSc::Clear() {
        local_size = global_size = 0;
        if (matrix != NULL) {
            MatrixDestroyDataPetsc(&matrix);
        }
        return true;
    }

    bool SolverPETSc::isMatrixSet() {
        return matrix != NULL;
    }

    std::string SolverPETSc::GetParameter(std::string name) const {
        std::cout << "SolverPETSc::GetParameter unsupported operation" << std::endl;
        //throw INMOST::SolverUnsupportedOperation;
        return "";
    }

    void SolverPETSc::SetParameter(std::string name, std::string value) {
        std::cout << "SolverPETSc::SetParameter unsupported operation" << std::endl;
        //throw INMOST::SolverUnsupportedOperation;
    }

    const INMOST_DATA_ENUM_TYPE SolverPETSc::Iterations() const {
        return static_cast<INMOST_DATA_ENUM_TYPE>(SolverIterationNumberPetsc(ksp));
    }

    const INMOST_DATA_REAL_TYPE SolverPETSc::Residual() const {
        return SolverResidualNormPetsc(ksp);
    }

    const std::string SolverPETSc::ReturnReason() const {
        return std::string(SolverConvergedReasonPetsc(ksp));
    }

    const std::string SolverPETSc::SolverName() const {
        return "petsc";
    }

    void SolverPETSc::Finalize() {
        if (petscSolversCount == 1) {
            SolverFinalizePetsc();
        }
    }

    SolverPETSc::~SolverPETSc() {
        petscSolversCount--;
    }

}
