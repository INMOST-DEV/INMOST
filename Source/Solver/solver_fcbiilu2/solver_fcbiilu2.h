#ifndef SOLVER_FCBIILU2_H_INCLUDED
#define SOLVER_FCBIILU2_H_INCLUDED
#include "inmost_solver.h"


//#define USE_SOLVER_FCBIILU2 //just for test; see flag HAVE_SOLVER_FCBIILU2
//#if defined(USE_SOLVER_FCBIILU2)
void MatrixCopyDataFcbiilu2(void ** ppA, void * pB);
void MatrixAssignDataFcbiilu2(void * pA, void * pB);
void MatrixInitDataFcbiilu2(void ** ppA, INMOST_MPI_Comm comm, const char * name);
void MatrixDestroyDataFcbiilu2(void ** pA);
void MatrixFillFcbiilu2(void * pA, int size, int nproc, int * ibl, int * ia, int * ja, double * values);
void MatrixFillValuesFcbiilu2(void * pA, double * values);
void MatrixFinalizeFcbiilu2(void * data);

void VectorInitDataFcbiilu2(void ** ppA, INMOST_MPI_Comm comm, const char * name);
void VectorCopyDataFcbiilu2(void ** ppA, void * pB);
void VectorAssignDataFcbiilu2(void * pA, void * pB);
void VectorPreallocateFcbiilu2(void * pA, int size);
void VectorFillFcbiilu2(void * pA, double * values);
void VectorLoadFcbiilu2(void * pA, double * values);
void VectorFinalizeFcbiilu2(void * data);
void VectorDestroyDataFcbiilu2(void ** ppA);

void SolverInitializeFcbiilu2(int * argc,char *** argv, const char * file_options);
bool SolverIsFinalizedFcbiilu2();
void SolverFinalizeFcbiilu2();
void SolverDestroyDataFcbiilu2(void ** data);
void SolverInitDataFcbiilu2(void ** data, INMOST_MPI_Comm comm, const char * name);
void SolverCopyDataFcbiilu2(void ** data, void * other_data, INMOST_MPI_Comm comm);
void SolverAssignDataFcbiilu2(void * data, void * other_data);
void SolverSetMatrixFcbiilu2(void * data, void * matrix_data, bool same_pattern, bool reuse_preconditioner);
bool SolverSolveFcbiilu2(void * data, void * rhs_data, void * sol_data);
int SolverIterationNumberFcbiilu2(void * data);
double SolverResidualNormFcbiilu2(void * data);
void SolverAddOtherStatFcbiilu2(void * data, unsigned int * pivmod, double * prdens, double * t_prec, double * t_iter);
//#endif //USE_SOLVER_FCBIILU2


#endif //SOLVER_FCBIILU2_H_INCLUDED
