#ifndef SOLVER_ANI_H_INCLUDED
#define SOLVER_ANI_H_INCLUDED
#include "inmost_solver.h"

void MatrixCopyDataAni(void ** ppA, void * pB);
void MatrixAssignDataAni(void * pA, void * pB);
void MatrixInitDataAni(void ** ppA, INMOST_MPI_Comm comm, const char * name);
void MatrixDestroyDataAni(void ** pA);
void MatrixFillAni(void * pA, int size, int * ia, int * ja, double * values);
void MatrixFillValuesAni(void * pA, double * values);
void MatrixFinalizeAni(void * data);

void VectorInitDataAni(void ** ppA, INMOST_MPI_Comm comm, const char * name);
void VectorCopyDataAni(void ** ppA, void * pB);
void VectorAssignDataAni(void * pA, void * pB);
void VectorPreallocateAni(void * pA, int size);
void VectorFillAni(void * pA, double * values);
void VectorLoadAni(void * pA, double * values);
void VectorFinalizeAni(void * data);
void VectorDestroyDataAni(void ** ppA);

void SolverInitializeAni(int * argc,char *** argv, const char * file_options);
bool SolverIsFinalizedAni();
void SolverFinalizeAni();
void SolverDestroyDataAni(void ** data);
void SolverInitDataAni(void ** data, INMOST_MPI_Comm comm, const char * name);
void SolverCopyDataAni(void ** data, void * other_data, INMOST_MPI_Comm comm);
void SolverAssignDataAni(void * data, void * other_data);
void SolverSetMatrixAni(void * data, void * matrix_data, bool same_pattern, bool reuse_preconditioner);
bool SolverSolveAni(void * data, void * rhs_data, void * sol_data);
int SolverIterationNumberAni(void * data);
double SolverResidualNormAni(void * data);


#endif //SOLVER_ANI_H_INCLUDED
