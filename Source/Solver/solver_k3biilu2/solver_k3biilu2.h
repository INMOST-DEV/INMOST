#ifndef SOLVER_K3BIILU2_H_INCLUDED
#define SOLVER_K3BIILU2_H_INCLUDED

#include "k3d_slv.hxx"
#include "inmost.h"
#include "inmost_solver.h"

typedef struct {
    int ittype;       // 0 - BiCGStab; 1,2,3 - GMRES(niter_cycle); 2 - +Poly(ncoef)_BiCGStab; 3 - +Poly(ncoef)_GMRESb(niter_cycle2)
    int niter_cycle;  // outer GMRES cycle (=kgmr by IEK); 1 - BiCGStab
    int ncoef;        // polynomial degree (=kdeg by IEK); ittype=2,3
    int niter_cycle2; // internal GMRES cycle size for ittype=3
    double eps;          // the residual precision: ||r|| < eps * ||b||; eps=1e-6
    int maxit;        // number of iterations permitted; maxit=999
    int ichk;         // number of skipped iterations to check the convergence
    int msglev;       // messages level; msglev=0 for silent; msglev>0 to output solution statistics
} ParIter;

void ParametersDefault(ParIter &parIter, k3d::SParams &parPrec);

/* BiCGStab solver structure */
typedef struct {
    int n;         // local number of unknowns at the local processor
    int nproc;     // total number of processors
    int *ibl;     // block splitting: ibl[0]=0; ibl[nproc]=nglob
    int *ia;      // row pointers: ia[0]=0; ia[nloc]=nzloc
    int *ja;      // global column numbers (NOTE: starting from 0 or 1); ja[nzloc]
    double *a;       // nonzero coefficients of local matrix A; a[nzloc]
    k3d::SParams *pParams;  // preconditioner construction parameters
    ParIter *pParIter; // preconditioner construction parameters
    k3d::CK3D_Solver<int, double, double> *pSolver;  // pointer to the solver structure
    int ierr;      // error flag on return; ierr=0 for success
    int istat[16]; // integer statistics array on return
    double dstat[16]; // double  statistics array on return
    double RESID;     // residual norm
    int ITER;      // number of BiCGStab iterations performed
    bool parameters_initialized; // backward compatibility with old parameter files
} bcg_k3biilu2;

typedef struct {
    int n;         // local number of unknowns at the local processor
    int nproc;     // total number of processors
    int *ibl;     // block splitting: ibl[0]=0; ibl[nproc]=nglob
    int *ia;      // row pointers: ia[0]=0; ia[nloc]=nzloc
    int *ja;      // global column numbers (NOTE: starting from 0 or 1); ja[nzloc]
    double *A;       // nonzero coefficients of local matrix A; a[nzloc]
} matrix_k3biilu2;

typedef struct {
    int n;         // local number of unknowns at the local processor
    double *v;       // local data vector
} vector_k3biilu2;

/*****************************************************************************/

/* Initialize bcg solver */
int initbcg_k3(bcg_k3biilu2 *s, matrix_k3biilu2 *A, double eps);

int newmatrixbcg_k3(bcg_k3biilu2 *s, matrix_k3biilu2 *A, bool same_precond);

/* Reinitialize solver preconditioner with new matrix A */
int renewbcg_k3(bcg_k3biilu2 *s, double *A);

/* Solve linear system */
int solvebcg_k3(bcg_k3biilu2 *s, vector_k3biilu2 *b, vector_k3biilu2 *x);

/* Free memory used by solver */
void freebcg_k3(bcg_k3biilu2 *s);

/*****************************************************************************/

//#define USE_SOLVER_K3BIILU2 //just for test; see flag HAVE_SOLVER_K3BIILU2
//#if defined(USE_SOLVER_K3BIILU2)
void MatrixCopyDataK3biilu2(matrix_k3biilu2 **ppA, matrix_k3biilu2 *pB);

void MatrixAssignDataK3biilu2(matrix_k3biilu2 *pA, matrix_k3biilu2 *pB);

void MatrixInitDataK3biilu2(matrix_k3biilu2 **ppA, INMOST_MPI_Comm comm, const char *name);

void MatrixDestroyDataK3biilu2(matrix_k3biilu2 **pA);

void MatrixFillK3biilu2(matrix_k3biilu2 *pA, int size, int nproc, int *ibl, int *ia, int *ja, double *values);

void MatrixFillValuesK3biilu2(matrix_k3biilu2 *pA, double *values);

void MatrixFinalizeK3biilu2(matrix_k3biilu2 *data);

void SolverInitializeK3biilu2(bcg_k3biilu2 *data, int *argc, char ***argv, const char *file_options);

bool SolverIsFinalizedK3biilu2();

void SolverFinalizeK3biilu2();

void SolverDestroyDataK3biilu2(bcg_k3biilu2 **data);

void SolverInitDataK3biilu2(bcg_k3biilu2 **data, INMOST_MPI_Comm comm, const char *name);

void SolverCopyDataK3biilu2(bcg_k3biilu2 **data, bcg_k3biilu2 *other_data, INMOST_MPI_Comm comm);

void SolverAssignDataK3biilu2(bcg_k3biilu2 *data, bcg_k3biilu2 *other_data);

void SolverSetMatrixK3biilu2(bcg_k3biilu2 *data, matrix_k3biilu2 *matrix_data, bool same_pattern, bool reuse_preconditioner);

bool SolverSolveK3biilu2(bcg_k3biilu2 *data, vector_k3biilu2 *rhs_data, vector_k3biilu2 *sol_data);


#endif //SOLVER_K3BIILU2_H_INCLUDED
