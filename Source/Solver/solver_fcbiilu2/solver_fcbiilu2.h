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
} bcg;

typedef struct
{
    int n;            // local number of unknowns at the local processor
    int nproc;        // number of processors
    int * ibl;        // block splitting: ibl[0]=0; ibl[nproc]=nglob
    int * ia;
    int * ja;
    double * A;
} matrix;

typedef struct
{
    int n;            // local number of unknowns at the local processor
    double * v;
} vector;

void MatrixCopyDataFcbiilu2(matrix **pA, matrix *B);
void MatrixAssignDataFcbiilu2(matrix *A, matrix* B);
void MatrixInitDataFcbiilu2(matrix ** ppA, INMOST_MPI_Comm comm, const char * name);
void MatrixDestroyDataFcbiilu2(matrix ** pA);
void MatrixFillFcbiilu2(matrix * pA, int size, int nproc, int * ibl, int * ia, int * ja, double * values);
void MatrixFillValuesFcbiilu2(matrix * pA, double * values);
void MatrixFinalizeFcbiilu2(matrix * data);

void VectorInitDataFcbiilu2(vector ** ppA, INMOST_MPI_Comm comm, const char * name);
void VectorCopyDataFcbiilu2(vector ** ppA, vector * pB);
void VectorAssignDataFcbiilu2(vector * pA, vector * pB);
void VectorPreallocateFcbiilu2(vector * pA, int size);
void VectorFillFcbiilu2(vector * pA, double * values);
void VectorLoadFcbiilu2(vector * pA, double * values);
void VectorFinalizeFcbiilu2(vector * data);
void VectorDestroyDataFcbiilu2(vector ** ppA);

void SolverInitializeFcbiilu2(bcg *data, int * argc,char *** argv, const char * file_options);
bool SolverIsFinalizedFcbiilu2();
void SolverFinalizeFcbiilu2();
void SolverDestroyDataFcbiilu2(bcg ** data);
void SolverInitDataFcbiilu2(bcg ** data, INMOST_MPI_Comm comm, const char * name);
void SolverCopyDataFcbiilu2(bcg **data, bcg *other_data, INMOST_MPI_Comm comm);
void SolverAssignDataFcbiilu2(bcg * data, bcg * other_data);
void SolverSetMatrixFcbiilu2(bcg * data, matrix * matrix_data, bool same_pattern, bool reuse_preconditioner);
bool SolverSolveFcbiilu2(bcg * data, vector * rhs_data, vector * sol_data);


#endif //SOLVER_FCBIILU2_H_INCLUDED
