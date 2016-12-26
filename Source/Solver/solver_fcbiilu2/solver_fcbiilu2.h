#ifndef SOLVER_FCBIILU2_H_INCLUDED
#define SOLVER_FCBIILU2_H_INCLUDED
#include "inmost_solver.h"

/* BiCGStab solver structure */
typedef struct
{
    int n;            // local number of unknowns at the local processor
    int nproc;        // number of processors
    int    * ibl;     // block splitting: ibl[0]=0; ibl[nproc]=nglob
    int    * ia;      // row pointers: ia[0]=0; ia[nloc]=nzloc
    int    * ja;      // column numbers (NOTE: starting from 0 or 1); ja[nzloc]
    double * a;       // matrix A coefficients; a[nzloc]
    int len_r8;       // size of the working memory, set len_r8=nbl for the first call
    double * W;       // poiter to the working memory W[len_r8]
    int    kovl;      // number of overlap layers: kovl=0,1,2,...
    double tau;       // the ILU2 precision (for the submatrix factorization); tau=3e-3
    double eps;       // the residual precision: ||r|| < eps * ||b||; eps=1e-6
    int    nit;       // number of iterations permitted; nit=999
    int    msglev;    // messages level; msglev=0 for silent; msglev=1 to output solution statistics
    int    ierr;      // error flag on return; ierr=0 for success
    int    istat[16]; // integer statistics array on return
    double dstat[16]; // double  statistics array on return
    double RESID;     // residual norm
    int    ITER;      // number of BiCGStab iterations performed
} bcg_fcbiilu2;

typedef struct
{
    int n;            // local number of unknowns at the local processor
    int nproc;        // number of processors
    int * ibl;        // block splitting: ibl[0]=0; ibl[nproc]=nglob
    int * ia;
    int * ja;
    double * A;
} matrix_fcbiilu2;

typedef struct
{
    int n;            // local number of unknowns at the local processor
    double * v;
} vector_fcbiilu2;

void MatrixCopyDataFcbiilu2(matrix_fcbiilu2 **pA, matrix_fcbiilu2 *B);
void MatrixAssignDataFcbiilu2(matrix_fcbiilu2 *A, matrix_fcbiilu2* B);
void MatrixInitDataFcbiilu2(matrix_fcbiilu2 ** ppA, INMOST_MPI_Comm comm, const char * name);
void MatrixDestroyDataFcbiilu2(matrix_fcbiilu2 ** pA);
void MatrixFillFcbiilu2(matrix_fcbiilu2 * pA, int size, int nproc, int * ibl, int * ia, int * ja, double * values);
void MatrixFillValuesFcbiilu2(matrix_fcbiilu2 * pA, double * values);
void MatrixFinalizeFcbiilu2(matrix_fcbiilu2 * data);

void VectorInitDataFcbiilu2(vector_fcbiilu2 ** ppA, INMOST_MPI_Comm comm, const char * name);
void VectorCopyDataFcbiilu2(vector_fcbiilu2 ** ppA, vector_fcbiilu2 * pB);
void VectorAssignDataFcbiilu2(vector_fcbiilu2 * pA, vector_fcbiilu2 * pB);
void VectorPreallocateFcbiilu2(vector_fcbiilu2 * pA, int size);
void VectorFillFcbiilu2(vector_fcbiilu2 * pA, double * values);
void VectorLoadFcbiilu2(vector_fcbiilu2 * pA, double * values);
void VectorFinalizeFcbiilu2(vector_fcbiilu2 * data);
void VectorDestroyDataFcbiilu2(vector_fcbiilu2 ** ppA);

void SolverInitializeFcbiilu2(bcg_fcbiilu2 *data, int * argc,char *** argv, const char * file_options);
bool SolverIsFinalizedFcbiilu2();
void SolverFinalizeFcbiilu2();
void SolverDestroyDataFcbiilu2(bcg_fcbiilu2 ** data);
void SolverInitDataFcbiilu2(bcg_fcbiilu2 ** data, INMOST_MPI_Comm comm, const char * name);
void SolverCopyDataFcbiilu2(bcg_fcbiilu2 **data, bcg_fcbiilu2 *other_data, INMOST_MPI_Comm comm);
void SolverAssignDataFcbiilu2(bcg_fcbiilu2 * data, bcg_fcbiilu2 * other_data);
void SolverSetMatrixFcbiilu2(bcg_fcbiilu2 * data, matrix_fcbiilu2 * matrix_data, bool same_pattern, bool reuse_preconditioner);
bool SolverSolveFcbiilu2(bcg_fcbiilu2 * data, vector_fcbiilu2 * rhs_data, vector_fcbiilu2 * sol_data);


#endif //SOLVER_FCBIILU2_H_INCLUDED
