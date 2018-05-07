#include "SolverPETSc.h"

namespace INMOST {

    unsigned int SolverPETSc::petscSolversCount = 0;

    SolverPETSc::SolverPETSc() {
        petscSolversCount++;
        this->ksp = NULL;
        this->matrix = NULL;
        maximum_iterations = 2500;
        schwartz_overlap = 1;

        atol = 1.0e-5;
        rtol = 1.0e-12;
        dtol = 1.0e+100;
        drop_tolerance = 0.005;
        fill_level = 3;
    }

    SolverInterface *SolverPETSc::Copy(const SolverInterface *other) {
        if (other == NULL) {
            throw INMOST::SolverCopyNullException;
        }
        const SolverPETSc *solver;
        try {
            solver = dynamic_cast<const SolverPETSc *>(other);
        } catch (...) {
            throw INMOST::SolverCopyException;
        }
        this->ksp = NULL;
        this->matrix = NULL;

        SolverCopyDataPetsc(&ksp, solver->ksp, solver->communicator);
        if (solver->matrix != NULL) {
            MatrixCopyDataPetsc(&matrix, solver->matrix);
            SolverSetMatrixPetsc(ksp, matrix, false, false);
        }
        return this;
    }

    void SolverPETSc::Assign(const SolverInterface *other) {
        if (other == NULL) {
            throw INMOST::SolverAssignNullException;
        }
        const SolverPETSc *other_solver;
        try {
            other_solver = dynamic_cast<const SolverPETSc *>(other);
        } catch (...) {
            throw INMOST::SolverAssignException;
        }
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

    void SolverPETSc::Setup(int *argc, char ***argv, SolverParameters &p) {
        this->parametersFile = p.internalFile;
        if (p.internalFile.empty()) {
            for (parameters_iterator_t parameter = p.parameters.begin(); parameter < p.parameters.end(); parameter++) {
                this->SetParameter((*parameter).first, (*parameter).second);
            }
        }

        SolverInitializePetsc(argc, argv, p.internalFile.empty() ? NULL : p.internalFile.c_str());
        SolverInitDataPetsc(&ksp, this->communicator, p.solverPrefix.c_str());
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
            SolverSetDropTolerancePetsc(ksp, drop_tolerance);
            SolverSetFillLevelPetsc(ksp, fill_level);
            SolverSetOverlapPetsc(ksp, schwartz_overlap);
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
            SolverSetTolerancesPetsc(ksp, rtol, atol, dtol, maximum_iterations);
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
	if (ksp != NULL) {
        	SolverDestroyDataPetsc(&ksp);
	}
        return true;
    }

    bool SolverPETSc::isMatrixSet() {
        return matrix != NULL;
    }

    std::string SolverPETSc::GetParameter(std::string name) const {
        if (name == "maximum_iterations") return to_string(maximum_iterations);
        else if (name == "schwartz_overlap") return to_string(schwartz_overlap);
        else if (name == "absolute_tolerance") return to_string(atol);
        else if (name == "relative_tolerance") return to_string(rtol);
        else if (name == "divergence_tolerance") return to_string(dtol);
        else if (name == "drop_tolerance") return to_string(drop_tolerance);
        else if (name == "fill_level") return to_string(fill_level);
        else {
#if !defined(SILENCE_SET_PARAMETER)
            std::cout << "Parameter " << name << " is unknown (Use internal file for all parameters)" << std::endl;
#endif
            return "";
        }
    }

    void SolverPETSc::SetParameter(std::string name, std::string value) {
        const char *val = value.c_str();
        if (name == "maximum_iterations") maximum_iterations = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "schwartz_overlap") schwartz_overlap = static_cast<INMOST_DATA_ENUM_TYPE>(atoi(val));
        else if (name == "absolute_tolerance") atol = atof(val);
        else if (name == "relative_tolerance") rtol = atof(val);
        else if (name == "divergence_tolerance") dtol = atof(val);
        else if (name == "drop_tolerance") drop_tolerance = atof(val);
        else if (name == "fill_level") fill_level = atof(val);
#if !defined(SILENCE_SET_PARAMETER)
        else std::cout << "Parameter " << name << " is unknown (Use internal file for all parameters)" << std::endl;
#endif
    }

    INMOST_DATA_ENUM_TYPE SolverPETSc::Iterations() const {
        return static_cast<INMOST_DATA_ENUM_TYPE>(SolverIterationNumberPetsc(ksp));
    }

    INMOST_DATA_REAL_TYPE SolverPETSc::Residual() const {
        return SolverResidualNormPetsc(ksp);
    }

    const std::string SolverPETSc::ReturnReason() const {
        return std::string(SolverConvergedReasonPetsc(ksp));
    }

    const std::string SolverPETSc::SolverName() const {
        return "petsc";
    }

    void SolverPETSc::Finalize() {
       // if (petscSolversCount == 1) {
        //    SolverFinalizePetsc(); \\commented by Kramarenko
        //}
    }

    SolverPETSc::~SolverPETSc() {
        this->Clear();
        //petscSolversCount--;
    }

}
