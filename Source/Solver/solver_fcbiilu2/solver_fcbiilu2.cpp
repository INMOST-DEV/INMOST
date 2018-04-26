#include "solver_fcbiilu2.h"

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <string>

#define T(x) // x // Trace of function calls. Use: "T(x) x" for trace and "T(x)" for silence

/*****************************************************************************/

//----- File: biilu2.h if available
// Solve the linear system A X = B by
//  biilu2_bcg solver with working memory allocations and statistics output
int biilu2_bcg(
        int *ibl,   // block splitting: ibl[0]=0; ibl[nproc]=n
        int *ia,    // row pointers: ia[0]=0; ia[nloc]=nzloc
        int *ja,    // column numbers (NOTE: starting from 0 or 1); ja[nzloc]
        double *a,  // matrix A coefficients; a[nzloc]
        double *b,  // right-hand side B; b[nloc]
        double *x,  // initial guess to the solution X on entry; solution X on return; x[nloc]
        int job,    // job number: 0 - construct preconditioner; 1 - use the previously constructed one
        int *len_r8,// size of the working memory, set len_r8=nbl for the first call
        double **S, // poiter to the working memory S[len_r8]
        int kovl,   // number of overlap layers: kovl=0,1,2,...
        double tau, // the ILU2 precision (for the submatrix factorization); tau=3e-3
        double eps, // the residual precision: ||r|| < eps * ||b||; eps=1e-6
        int nit,    // number of iterations permitted; nit=999; if(nit==0) preconditioner construction only
        int msglev, // messages level; msglev=0 for silent; msglev=1 to output solution statistics
        int *ierr,  // error flag on return; ierr=0 for success
        int *istat,  //[16],  // integer statistics array on return
        double *dstat); //[16]); // double  statistics array on return
//
// Here, in notation:
//      nproc - the number of processors (or blocks), nproc is equal to the MPI communicator size
//      nzloc - the local number of nonzero elements, nzloc=ia[[myid+1]]-ia[ibl[myid]]
//      nloc  - the local number of unknowns at the current processor, nloc=ibl[myid+1]-ibl[myid]
//      n     - the total number of unknowns in matrix A, n=ibl[nproc]
//
// ON ENTRY:
//      ibl, ia, ja, a, b, x, job, len_r8, S, kovl, tau, eps, nit, msglev
// ON RETURN:
//      x, len_r8, S, ierr, istat, dstat
//
// Run selftest by
//      mpicc -lm -DBIILU2_BICGST_SELFTEST biilu2.c && mpirun -np 2 a.out

/*****************************************************************************/

/* Initialize bcg solver */
static int initbcg(bcg_fcbiilu2 *s, matrix_fcbiilu2 *A, double eps);

/* Reinitialize solver preconditioner with new matrix A */
static int renewbcg(bcg_fcbiilu2 *s, double *A);
/* Solve linear system */
/*static*/ int solvebcg(bcg_fcbiilu2 *s, vector_fcbiilu2 *b, vector_fcbiilu2 *x);

/* Free memory used by solver */
static void freebcg(bcg_fcbiilu2 *s);

/*****************************************************************************/

/* Initialize solver with new matrix A */
static int newmatrixbcg(bcg_fcbiilu2 *s, matrix_fcbiilu2 *A, bool same_precond) {
    if (s->n != A->n && same_precond) throw INMOST::CannotReusePreconditionerOfDifferentSize;
    s->n = A->n;
    s->nproc = A->nproc;
    s->ibl = A->ibl;
    s->ia = A->ia;
    s->ja = A->ja;
    if (!same_precond) {
        //do nothing...
        T(std::cout << "##### inside newmatrixbcg bef. renewbcg \n";)//db!
        return renewbcg(s, A->A);
    } else return 0;
}

/* solver */
/* Initialize bcg solver */
int initbcg(bcg_fcbiilu2 *s, matrix_fcbiilu2 *A) {
    //s->eps = eps;
    T(std::cout << "##### inside initbcg bef. newmatrixbcg eps=" << eps << " \n";)//db!
    return newmatrixbcg(s, A, false);
}

/* Reinitialize solver preconditioner with new matrix A */
int renewbcg(bcg_fcbiilu2 *s, double *A) {
    //reinitialize matrix values
    s->a = A;
    //BIILU2 preconditioner construction...
    //s->kovl = set_kovl;
    //s->tau = set_tau;
    //s->msglev = set_msglev;

    int job = 0;
    int maxit = 0;
    s->len_r8 = 0; //to be set by nbl inside biilu2_bcg
    if (s->W) free(s->W);
    s->W = NULL;

    int ierr = 0;
    T(std::cout << "##### inside renewbcg bef. biilu2_bcg\n";)//db!
    double *B = (double *) malloc(sizeof(double) * s->n); //db!!!!!!!!!!!!!!
    //double *X = (double*) malloc(sizeof(double)*s->n); //db!!!!!!!!!!!!!!
    biilu2_bcg(s->ibl, s->ia, s->ja, s->a,
            //B, X, //!!!!!KAPORIN!!!!!!!!!!!!!
               B, NULL, //!!!!!KAPORIN!!!!!!!!!!!!
            //NULL, NULL, //!!!!!INK!!!!!!!!!!!
               job, &s->len_r8, &s->W,
               s->kovl, s->tau, s->eps, maxit, s->msglev,
               &ierr, s->istat, s->dstat);
    T(std::cout << "##### inside renewbcg aft. biilu2_bcg\n";)//db!
    free(B); //free(X);//db!!!!!!!!!!!!!!!
    if (ierr) printf("initialization of biilu2 failed, ierr=%d\n", ierr);

    return ierr;
}

/* Solve linear system */
int solvebcg(bcg_fcbiilu2 *s, vector_fcbiilu2 *b, vector_fcbiilu2 *x) {
    //s->kovl = set_kovl;
    //s->tau = set_tau;
//  s->eps     = set_eps;
    //s->nit = set_nit;
    //s->msglev = set_msglev;

    int job = 1;
    int maxit = s->nit;

    int ierr = 0;
    T(std::cout << "##### inside solvebcg bef. biilu2_bcg\n";)//db!
    biilu2_bcg(s->ibl, s->ia, s->ja, s->a, b->v, x->v,
               job, &s->len_r8, &s->W,
               s->kovl, s->tau, s->eps, maxit, s->msglev,
               &ierr, s->istat, s->dstat);
    T(std::cout << "##### inside solvebcg aft. biilu2_bcg\n";)//db!

    s->ITER = s->istat[2];
    s->RESID = s->dstat[2];

    return ierr;
}

/* Free memory used by solver */
void freebcg(bcg_fcbiilu2 *s) {
    if (s->W) free(s->W);
    s->W = NULL;
}

/*****************************************************************************/

void MatrixCopyDataFcbiilu2(matrix_fcbiilu2 **pA, matrix_fcbiilu2 *B) {
    if (pA == NULL || B == NULL) throw INMOST::DataCorruptedInSolver;
    *pA = (matrix_fcbiilu2 *) malloc(sizeof(matrix_fcbiilu2));
    matrix_fcbiilu2 *A = *pA;
    A->n = B->n;
    if (B->n != 0) {
        int nnz = B->ia[B->n] - B->ia[0];
        A->nproc = B->nproc;
        A->ibl = (int *) malloc(sizeof(int) * (A->nproc + 1));
        memcpy(A->ibl, B->ibl, sizeof(int) * (A->nproc + 1));
        A->ia = (int *) malloc(sizeof(int) * (A->n + 1));
        memcpy(A->ia, B->ia, sizeof(int) * (A->n + 1));
        A->ja = (int *) malloc(sizeof(int) * nnz);
        memcpy(A->ja, B->ja, sizeof(int) * nnz);
        A->A = (double *) malloc(sizeof(double) * nnz);
        memcpy(A->A, B->A, sizeof(double) * nnz);
    }
}

void MatrixAssignDataFcbiilu2(matrix_fcbiilu2 *A, matrix_fcbiilu2 *B) {
    if (A == NULL || B == NULL) throw INMOST::DataCorruptedInSolver;
    if (A != B) {
        if (A->n != 0) {
            free(A->ibl);
            free(A->ia);
            free(A->ja);
            free(A->A);
        }
        if (B->n != 0) {
            int nnz = B->ia[B->n] - B->ia[0];
            A->n = B->n;
            A->nproc = B->nproc;
            A->ibl = (int *) malloc(sizeof(int) * (A->nproc + 1));
            memcpy(A->ibl, B->ibl, sizeof(int) * (A->nproc + 1));
            A->ia = (int *) malloc(sizeof(int) * (A->n + 1));
            memcpy(A->ia, B->ia, sizeof(int) * (A->n + 1));
            A->ja = (int *) malloc(sizeof(int) * nnz);
            memcpy(A->ja, B->ja, sizeof(int) * nnz);
            A->A = (double *) malloc(sizeof(double) * nnz);
            memcpy(A->A, B->A, sizeof(double) * nnz);
        }
    }
}

void MatrixInitDataFcbiilu2(matrix_fcbiilu2 **ppA, INMOST_MPI_Comm comm, const char *name) {
    T(std::cout << "##### ins. MatrixInitDataFcbiilu2 \n";)//db!
    if (ppA == NULL) throw INMOST::DataCorruptedInSolver;
    if (*ppA == NULL) {
        *ppA = (matrix_fcbiilu2 *) malloc(sizeof(matrix_fcbiilu2));
        matrix_fcbiilu2 *A = (matrix_fcbiilu2 *) *ppA;
        A->n = 0;
        A->nproc = 0;
        T(std::cout << "##### ins. MatrixInitDataFcbiilu2 n=nproc=0 \n";)//db!
    }
    (void) comm;
    (void) name;
}

void MatrixDestroyDataFcbiilu2(matrix_fcbiilu2 **pA) {
    matrix_fcbiilu2 *A = (matrix_fcbiilu2 *) (*pA);
    if (A != NULL) {
        if (A->n != 0) {
            free(A->ibl);
            free(A->ia);
            free(A->ja);
            free(A->A);
            T(std::cout << "##### ins. MatrixDestroyDataFcbiilu2 ...free \n";)//db!
        }
        free(*pA);
        *pA = NULL;
    }
}


void MatrixFillFcbiilu2(matrix_fcbiilu2 *pA, int size, int nproc, int *ibl, int *ia, int *ja, double *values) {
    T(std::cout << "##### ins. MatrixFillFcbiilu2 n=" << size << " nproc=" << nproc << " \n";)//db!
    if (pA == NULL) throw INMOST::DataCorruptedInSolver;
    matrix_fcbiilu2 *A = (matrix_fcbiilu2 *) pA;
    A->n = size;
    A->nproc = nproc;
    A->ibl = ibl;
    A->ia = ia;
    A->ja = ja;
    A->A = values;
}

void MatrixFillValuesFcbiilu2(matrix_fcbiilu2 *pA, double *values) {
    T(std::cout << "##### ins. MatrixFillValuesFcbiilu2 \n";)//db!
    if (pA == NULL) throw INMOST::DataCorruptedInSolver;
    matrix_fcbiilu2 *A = (matrix_fcbiilu2 *) pA;
    free(A->A);
    A->A = values;
}

void MatrixFinalizeFcbiilu2(matrix_fcbiilu2 *data) {
    //don't need to do anything
    (void) data;
}

void SolverInitializeFcbiilu2(bcg_fcbiilu2 *data, int *argc, char ***argv, const char *file_options) {
    T(std::cout << "##### ins. SolverInitializeFcbiilu2 (" << file_options << ") \n";)//db!
    if (file_options == NULL) return;
    std::string s = file_options;
    if (s == "" || s == " ") return;
    std::ifstream is;
    T(std::cout << "##### ins. SolverInitializeFcbiilu2: bef. open(" << file_options << ") \n";)//db!
    is.open(file_options, std::ifstream::in);
    if (s == "ctrl_dat") {
        getline(is, s);                                      //1 skip ipart
        getline(is, s);                                      //2 skip mtx filename
        getline(is, s);                                      //3 skip rhs filename
        getline(is, s);
        sscanf(s.c_str(), "%lg", &(data->eps));  //4 eps
        getline(is, s);
        sscanf(s.c_str(), "%d", &(data->nit));  //5 nit
        getline(is, s); //sscanf(s.c_str(), "%d",  &set_kgmr); //6- skip kgmr
        getline(is, s); //sscanf(s.c_str(), "%d",  &set_kdeg); //7- skip kdeg
        getline(is, s);
        sscanf(s.c_str(), "%d", &(data->kovl)); //8 kovl
        getline(is, s);
        sscanf(s.c_str(), "%lg", &(data->tau));  //9 tau
        //? msglev
        data->params_initialized = true;
    } else if (s == "biilu2_options.txt") { // file: "biilu2_options.txt"
        getline(is, s);
        sscanf(s.c_str(), "%d", &(data->kovl));   //1 kovl
        getline(is, s);
        sscanf(s.c_str(), "%lg", &(data->tau));    //2 tau
        getline(is, s);
        sscanf(s.c_str(), "%lg", &(data->eps));    //3 eps
        getline(is, s);
        sscanf(s.c_str(), "%d", &(data->nit));    //4 nit
        getline(is, s);
        sscanf(s.c_str(), "%d", &(data->msglev)); //5 msglev
        data->params_initialized = true;
    }
    T(std::cout << "##### ins. SolverInitializeFcbiilu2:  kovl=" << set_kovl << " tau=" << set_tau << " eps=" << set_eps << " nit=" << set_nit << " msglev="
                << set_msglev << " from: " << file_options << " \n";)//db!
    (void) argc;
    (void) argv;
}

bool SolverIsFinalizedFcbiilu2() {
    return true; //no need to finalize
}

void SolverFinalizeFcbiilu2() {
}

void SolverDestroyDataFcbiilu2(bcg_fcbiilu2 **data) {
    if (data != NULL) {
        if (*data != NULL) {
            bcg_fcbiilu2 *m = (bcg_fcbiilu2 *) *data;
            freebcg(m);
            free(m);
        }
        *data = NULL;
    }
}

void SolverInitDataFcbiilu2(bcg_fcbiilu2 **data, INMOST_MPI_Comm comm, const char *name) {
    T(std::cout << "##### ins. SolverInitDataFcbiilu2 \n";)//db!
    *data = (bcg_fcbiilu2 *) malloc(sizeof(bcg_fcbiilu2));
    ((bcg_fcbiilu2 *) *data)->n = 0;
    ((bcg_fcbiilu2 *) *data)->nproc = 0;
    ((bcg_fcbiilu2 *) *data)->len_r8 = 0;
    ((bcg_fcbiilu2 *) *data)->W = NULL;
    (void) comm;
    (void) name;
}

void SolverCopyDataFcbiilu2(bcg_fcbiilu2 **data, bcg_fcbiilu2 *other_data, INMOST_MPI_Comm comm) {
    throw INMOST::NotImplemented; //later
    (void) data;
    (void) other_data;
    (void) comm;
}

void SolverAssignDataFcbiilu2(bcg_fcbiilu2 *data, bcg_fcbiilu2 *other_data) {
    throw INMOST::NotImplemented; //later
    (void) data;
    (void) other_data;
}

void SolverSetMatrixFcbiilu2(bcg_fcbiilu2 *data, matrix_fcbiilu2 *matrix_data, bool same_pattern, bool reuse_preconditioner) {
    T(std::cout << "##### ins. SolverSetMatrixFcbiilu2 \n";)//db!
    bcg_fcbiilu2 *m = (bcg_fcbiilu2 *) data;
    matrix_fcbiilu2 *A = (matrix_fcbiilu2 *) matrix_data;
    T(if (A == NULL) std::cout << "##### A == NULL ... \n";)//db!
    T(if (m == NULL) std::cout << "##### m == NULL ... \n";)//db!
    if (A == NULL || m == NULL) throw INMOST::DataCorruptedInSolver;
    T(std::cout << "##### ins. SolverSetMatrixFcbiilu2 bef. initbcg or newmatrixbcg \n";)//db!
    if (m->n == 0)
        initbcg(m, A);
    else
        newmatrixbcg(m, A, reuse_preconditioner);
    (void) same_pattern;
    T(std::cout << "##### ins. SolverSetMatrixFcbiilu2 bef. return \n";)//db!
}

bool SolverSolveFcbiilu2(bcg_fcbiilu2 *data, vector_fcbiilu2 *rhs_data, vector_fcbiilu2 *sol_data) {
    T(std::cout << "##### ins. SolverSolveFcbiilu2 \n";)//db!
    bcg_fcbiilu2 *m = (bcg_fcbiilu2 *) data;
    vector_fcbiilu2 *rhs = (vector_fcbiilu2 *) rhs_data, *sol = (vector_fcbiilu2 *) sol_data;
    return solvebcg(m, rhs, sol) == 0;
}
