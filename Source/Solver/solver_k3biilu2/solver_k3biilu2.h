#ifndef SOLVER_K3BIILU2_H_INCLUDED
#define SOLVER_K3BIILU2_H_INCLUDED
#include "inmost_solver.h"


//#define USE_SOLVER_K3BIILU2 //just for test; see flag HAVE_SOLVER_K3BIILU2
//#if defined(USE_SOLVER_K3BIILU2)
void MatrixCopyDataK3biilu2(void ** ppA, void * pB);
void MatrixAssignDataK3biilu2(void * pA, void * pB);
void MatrixInitDataK3biilu2(void ** ppA, INMOST_MPI_Comm comm, const char * name);
void MatrixDestroyDataK3biilu2(void ** pA);
void MatrixFillK3biilu2(void * pA, int size, int nproc, int * ibl, int * ia, int * ja, double * values);
void MatrixFillValuesK3biilu2(void * pA, double * values);
void MatrixFinalizeK3biilu2(void * data);

void VectorInitDataK3biilu2(void ** ppA, INMOST_MPI_Comm comm, const char * name);
void VectorCopyDataK3biilu2(void ** ppA, void * pB);
void VectorAssignDataK3biilu2(void * pA, void * pB);
void VectorPreallocateK3biilu2(void * pA, int size);
void VectorFillK3biilu2(void * pA, double * values);
void VectorLoadK3biilu2(void * pA, double * values);
void VectorFinalizeK3biilu2(void * data);
void VectorDestroyDataK3biilu2(void ** ppA);

void SolverInitializeK3biilu2(int * argc,char *** argv, const char * file_options);
bool SolverIsFinalizedK3biilu2();
void SolverFinalizeK3biilu2();
void SolverDestroyDataK3biilu2(void ** data);
void SolverInitDataK3biilu2(void ** data, INMOST_MPI_Comm comm, const char * name);
void SolverCopyDataK3biilu2(void ** data, void * other_data, INMOST_MPI_Comm comm);
void SolverAssignDataK3biilu2(void * data, void * other_data);
void SolverSetMatrixK3biilu2(void * data, void * matrix_data, bool same_pattern, bool reuse_preconditioner);
bool SolverSolveK3biilu2(void * data, void * rhs_data, void * sol_data);
int SolverIterationNumberK3biilu2(void * data);
double SolverResidualNormK3biilu2(void * data);
void SolverAddOtherStatK3biilu2(void * data, unsigned int * pivmod, double * prdens, double * t_prec, double * t_iter);
//#endif //USE_K3SOLVER_BIILU2


#endif //SOLVER_K3BIILU2_H_INCLUDED
