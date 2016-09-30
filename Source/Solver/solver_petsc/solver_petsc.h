//
// Created by Dmitri Bagaev on 28.09.16.
//

#ifndef INMOST_NEW_SOLVER_PETSC_H
#define INMOST_NEW_SOLVER_PETSC_H

//IF CERTAIN FUNCTIONALITY IS NOT AVAILABLE, PLEASE THROW INMOST::NotImplemented EXCEPTION

#include <inmost_common.h>

void MatrixInitDataPetsc(Mat **matrix, INMOST_MPI_Comm comm, const char * name);
void MatrixCopyDataPetsc(Mat **matrix, Mat *other_matrix);
void MatrixAssignDataPetsc(Mat *matrix, Mat *other_matrix);
void MatrixDestroyDataPetsc(Mat **matrix);
void MatrixPreallocatePetsc(Mat *matrix, int local_size, int global_size, int * diag_nonzeroes, int * off_diag_nonzeroes);
void MatrixFillPetsc(Mat *matrix, int row, int cols, int * col_positions, double * col_values);
void MatrixFinalizePetsc(Mat *matrix);


void VectorInitDataPetsc(Vec **vector, INMOST_MPI_Comm comm, const char * name);
void VectorCopyDataPetsc(Vec **vector, Vec *other_vector);
void VectorAssignDataPetsc(Vec *vector, Vec *other_vector);
void VectorDestroyDataPetsc(Vec **vector);
void VectorPreallocatePetsc(Vec *vector, int local_size, int global_size);
void VectorFillPetsc(Vec *vector, int size, int * positions, double * values);
void VectorLoadPetsc(Vec *vector, int size, int * positions, double * values);
void VectorFinalizePetsc(Vec *vector);

bool SolverIsInitializedPetsc();
void SolverInitializePetsc(int * argc,char *** argv, const char * file_options);
bool SolverIsFinalizedPetsc();
void SolverFinalizePetsc();

void SolverInitDataPetsc(KSP **ksp, INMOST_MPI_Comm comm, const char * name);
void SolverCopyDataPetsc(KSP **ksp, KSP *other_ksp, INMOST_MPI_Comm comm);
void SolverAssignDataPetsc(KSP *ksp, KSP *other_ksp);
void SolverDestroyDataPetsc(KSP **ksp);
void SolverSetMatrixPetsc(KSP *ksp, Mat *matrix, bool same_pattern, bool reuse_preconditioner);
bool SolverSolvePetsc(KSP *ksp, Vec *rhs, Vec *sol);
int SolverIterationNumberPetsc(KSP *ksp);
double SolverResidualNormPetsc(KSP *ksp);
const char * SolverConvergedReasonPetsc(KSP *ksp);

void SolverSetTolerancesPetsc(KSP *ksp, double rtol, double atol, double divtol, int maxits);
void SolverSetOverlapPetsc(KSP *ksp, int levels);
void SolverSetDropTolerancePetsc(KSP *ksp, double dtol);
void SolverSetFillLevelPetsc(KSP *ksp, double lfill);

#endif //INMOST_NEW_SOLVER_PETSC_H
