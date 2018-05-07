#include "solver_k3biilu2.h"

//#if defined(USE_SOLVER_K3BIILU2)

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <string>

#include "k3d.h"

#define T(x) //x // Trace of function calls. Use: "T(x) x" for trace and "T(x)" for silence

//TODO: ncycle niter_cycle niter_cycle2 ncoef nit->maxit

/*****************************************************************************/
#include <stdlib.h>   // for malloc()

#if defined (USE_MPI)

#include <mpi.h>      // for MPI_COMM_WORLD etc.

#endif

void ParametersDefault(ParIter &parIter, k3d::SParams &parPrec) {
    parIter.ittype = 0;    // 0 - BiCGStab; 1,2,3 - GMRES(niter_cycle); 2 - +Poly(ncoef)_BiCGStab; 3 - +Poly(ncoef)_GMRESb(niter_cycle2)
    parIter.niter_cycle = 1;   // outer GMRES cycle (=kgmr by IEK); 1 - BiCGStab
    parIter.ncoef = 1;    // polynomial degree (=kdeg by IEK); ittype=2,3
    parIter.niter_cycle2 = 4;    // internal GMRES cycle size for ittype=3
    parIter.eps = 1e-6; // the residual precision: ||r|| < eps * ||b||; eps=1e-6
    parIter.maxit = 999;  // number of iterations permitted; maxit=999
    parIter.ichk = 5;    // number of skipped iterations to check the convergence
    parIter.msglev = 0;    // messages level; msglev=0 for silent; msglev>0 to output solution statistics

    parPrec.ncycle = 3;
    parPrec.tau1 = 3e-3;
    parPrec.tau2 = -1.0;
}

// Solve the linear system A X = B by
//  k3biilu2_bcg solver with working memory allocations and statistics output
int k3biilu2_bcg(
        int *ibl,   // block splitting: ibl[0]=0; ibl[nproc]=n
        int *ia,    // row pointers: ia[0]=0; ia[nloc]=nzloc
        int *ja,    // column numbers (NOTE: starting from 0 or 1); ja[nzloc]
        double *a,     // matrix A coefficients; a[nzloc]
        double *b,     // right-hand side B; b[nloc]
        double *x,     // initial guess to the solution X on entry; solution X on return; x[nloc]
        int job,    // job number: 0 - construct preconditioner; 1 - use the previously constructed one
        int maxit,  // number of iterations permitted; maxit=999; if(maxit==0) preconditioner construction only; it is more important than one in pParIter
        k3d::SParams *pParams,  // preconditioner construction parameters
        ParIter *pParIter, // iterative solver parameters
        k3d::CK3D_Solver<int, double, double> *pSolver,  // pointer to the solver structure
        int *ierr,  // error flag on return; ierr=0 for success
        int *istat, //[16],  // integer statistics array on return
        double *dstat) //[16]); // double  statistics array on return
//
// Here, in notation:
//      nproc - the number of processors (or blocks), nproc is equal to the MPI communicator size
//      nzloc - the local number of nonzero elements, nzloc=ia[[myid+1]]-ia[ibl[myid]]
//      nloc  - the local number of unknowns at the current processor, nloc=ibl[myid+1]-ibl[myid]
//      n     - the total number of unknowns in matrix A, n=ibl[nproc]
//
// ON ENTRY:
//      ibl, ia, ja, a, b, x, job, maxit, pParams, pParIter, pSolver
// ON RETURN:
//      x, pSolver, ierr, istat, dstat
//
{
    // Initialize MPI variables
    INMOST_MPI_Comm comm = INMOST_MPI_COMM_WORLD;
    int np = 1, mp = 0;
#if defined (USE_MPI)
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &mp);
#endif

    std::vector<long long> blks(np + 1);
    std::vector<int> blk2cpu(np + 1);

    long long *pblks = &blks[0];
    int *pblk2cpu = &blk2cpu[0];

    for (int i = 0; i <= np; i++) pblks[i] = ibl[i];
    for (int i = 0; i < np; i++) pblk2cpu[i] = i;

    T(cout << "HHHHH pSolver = " << pSolver << "\n";)//DB!

    if (mp != 0) pParIter->msglev = 0;

    if (job == 0) {
    } else {
        T(cout << "HHHHH k3biilu2_bcg: job == 0\n";)//DB!
        if (pParams->tau2 < -0.01) pParams->tau2 = pParams->tau1 * pParams->tau1;
        int nmodif;
        double prec_extend, density, scpiv_min, scpiv_max, piv_min, piv_max, dtime_fct;
        bool b_store_matrix = true;
        //if (maxit > 0) b_store_matrix = true;
        //TODO PIVMIN??? -> AUX = 0  <vars> / AUX
        //pParams->pivmin=0.0;
        pSolver->PrepareSolver((void *) &comm, pParams, np, pblks, pblk2cpu,
                               b_store_matrix, ia, ja, a,
                               prec_extend, density, scpiv_min, scpiv_max, nmodif, piv_min, piv_max, dtime_fct);

        if (pParIter->msglev > 0)
            std::cout << " K3: prec_extend=" << prec_extend << "  density=" << density << "  scpiv_min=" << scpiv_min << "  scpiv_max=" << scpiv_max
                      << "  nmodif=" << nmodif << "  piv_min=" << piv_min << "  piv_max=" << piv_max << "  dtime_fct=" << dtime_fct << endl;

        istat[0] = nmodif;
        dstat[0] = density;
        dstat[7] = dtime_fct;
//    } else {
//        T(cout<<"HHHHH k3biilu2_bcg: else (job == 0)\n";)//DB!
        //pSolver->PrepareMatrix((void *)&comm, np, pblks, pblk2cpu,
        //                       ia, ja, a);
    }

    if (maxit > 0) {
        T(cout << "HHHHH k3biilu2_bcg: 0<maxit=" << maxit << "\n";)//DB!
        ofstream *pfout = NULL;
        int niter, nmvm;
        double rhs_norm, res_ini, res_fin, dtime_iter;

        pSolver->SolveIter(pParIter->ittype, maxit, pParIter->niter_cycle, pParIter->ncoef, pParIter->niter_cycle2,
                           pParIter->eps, pParIter->ichk, pParIter->msglev, pfout,
                           b, x,
                           rhs_norm, res_ini, niter, nmvm, res_fin, dtime_iter);

        if (pParIter->msglev > 0)
            std::cout << " K3: rhs_norm=" << rhs_norm << "  res_ini=" << res_ini << "  niter=" << niter << "  nmvm=" << nmvm << "  res_fin=" << res_fin
                      << "  dtime_iter=" << dtime_iter << endl;

        dstat[2] = (rhs_norm == 0e0) ? res_fin : res_fin / rhs_norm;
        dstat[9] = dtime_iter;

        istat[2] = niter;
        pSolver->CleanMvmA(); //? prec iter iter
    }

    *ierr = 0; //TODO correct error status of construction or convergence

    return *ierr;
}

/* Initialize solver with new matrix A */
int newmatrixbcg_k3(bcg_k3biilu2 *s, matrix_k3biilu2 *A, bool same_precond) {
    if (s->n != A->n && same_precond) throw INMOST::CannotReusePreconditionerOfDifferentSize;
    s->n = A->n;
    s->nproc = A->nproc;
    s->ibl = A->ibl;
    s->ia = A->ia;
    s->ja = A->ja;
    if (!same_precond) {
        //do nothing...
        T(std::cout << "##### inside newmatrixbcg bef. renewbcg \n";)//db!
        return renewbcg_k3(s, A->A);
    } else return 0;
}

/* solver */
/* Initialize bcg solver */
int initbcg_k3(bcg_k3biilu2 *s, matrix_k3biilu2 *A, double eps) {
    s->pParIter->eps = eps;
    T(std::cout << "##### inside initbcg bef. newmatrixbcg eps=" << eps << " \n";)//db!
    return newmatrixbcg_k3(s, A, false);
}

/* Reinitialize solver preconditioner with new matrix A */
int renewbcg_k3(bcg_k3biilu2 *s, double *A) {
    //reinitialize matrix values
    s->a = A;
    //BIILU2 preconditioner construction...

    int job = 0;
    int maxit = 0;
    T(cout << "HHHHH bef. Clean in renewbcg\n";)//DB!
    s->pSolver->Clean();

    int ierr = 0;
    T(std::cout << "##### inside renewbcg bef. k3biilu2_bcg\n";)//db!
    k3biilu2_bcg(s->ibl, s->ia, s->ja, s->a,
                 NULL, NULL,
                 job, maxit, s->pParams, s->pParIter, s->pSolver,
                 &ierr, s->istat, s->dstat);
    T(std::cout << "##### inside renewbcg aft. k3biilu2_bcg ierr=" << ierr << "\n";)//db!
    if (ierr) std::cout << "initialization of k3biilu2 failed, ierr=" << ierr << endl; //TODO: myid

    return ierr;
}

/* Solve linear system */
int solvebcg_k3(bcg_k3biilu2 *s, vector_k3biilu2 *b, vector_k3biilu2 *x) {
    int job = 1;
    int maxit = s->pParIter->maxit;

    int ierr = 0;
    T(std::cout << "##### inside solvebcg bef. k3biilu2_bcg\n";)//db!
    k3biilu2_bcg(s->ibl, s->ia, s->ja, s->a, b->v, x->v,
                 job, maxit, s->pParams, s->pParIter, s->pSolver,
                 &ierr, s->istat, s->dstat);
    T(std::cout << "##### inside solvebcg aft. k3biilu2_bcg\n";)//db!
    if (ierr) std::cout << "linear system solution by k3biilu2 failed, ierr=" << ierr << endl; //TODO: myid

    s->ITER = s->istat[2];
    s->RESID = s->dstat[2];

    return ierr;
}

/* Free memory used by solver */
void freebcg_k3(bcg_k3biilu2 *s) {
    T(std::cout << "##### inside freebcg bef. clean\n";)//db!
    s->pSolver->Clean();
    delete s->pSolver;
    delete s->pParams;
    delete s->pParIter;
}

/*****************************************************************************/

void MatrixCopyDataK3biilu2(matrix_k3biilu2 **ppA, matrix_k3biilu2 *pB) {
    matrix_k3biilu2 *B = pB;
    if (ppA == NULL || pB == NULL) throw INMOST::DataCorruptedInSolver;
    *ppA = (matrix_k3biilu2 *) malloc(sizeof(matrix_k3biilu2));
    matrix_k3biilu2 *A = *ppA;
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

void MatrixAssignDataK3biilu2(matrix_k3biilu2 *pA, matrix_k3biilu2 *pB) {
    matrix_k3biilu2 *A = pA;
    matrix_k3biilu2 *B = pB;
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

void MatrixInitDataK3biilu2(matrix_k3biilu2 **ppA, INMOST_MPI_Comm comm, const char *name) {
    T(std::cout << "##### ins. MatrixInitDataK3biilu2 \n";)//db!
    if (ppA == NULL) throw INMOST::DataCorruptedInSolver;
    if (*ppA == NULL) {
        *ppA = (matrix_k3biilu2 *) malloc(sizeof(matrix_k3biilu2));
        matrix_k3biilu2 *A = *ppA;
        A->n = 0;
        A->nproc = 0;
        T(std::cout << "##### ins. MatrixInitDataK3biilu2 n=nproc=0 \n";)//db!
    }
    (void) comm;
    (void) name;
}

void MatrixDestroyDataK3biilu2(matrix_k3biilu2 **pA) {
    matrix_k3biilu2 *A = (*pA);
    if (A != NULL) {
        if (A->n != 0) {
            free(A->ibl);
            free(A->ia);
            free(A->ja);
            free(A->A);
            T(std::cout << "##### ins. MatrixDestroyDataK3biilu2 ...free \n";)//db!
        }
        free(*pA);
        *pA = NULL;
    }
}


void MatrixFillK3biilu2(matrix_k3biilu2 *pA, int size, int nproc, int *ibl, int *ia, int *ja, double *values) {
    T(std::cout << "##### ins. MatrixFillK3biilu2 n=" << size << " nproc=" << nproc << " \n";)//db!
    if (pA == NULL) throw INMOST::DataCorruptedInSolver;
    matrix_k3biilu2 *A = pA;
    A->n = size;
    A->nproc = nproc;
    A->ibl = ibl;
    A->ia = ia;
    A->ja = ja;
    A->A = values;
}

void MatrixFillValuesK3biilu2(matrix_k3biilu2 *pA, double *values) {
    T(std::cout << "##### ins. MatrixFillValuesK3biilu2 \n";)//db!
    if (pA == NULL) throw INMOST::DataCorruptedInSolver;
    matrix_k3biilu2 *A = pA;
    free(A->A);
    A->A = values;
}

void MatrixFinalizeK3biilu2(matrix_k3biilu2 *data) {
    //don't need to do anything
    (void) data;
}

void SolverInitializeK3biilu2(bcg_k3biilu2 *data, int *argc, char ***argv, const char *file_options) {
    T(std::cout << "##### ins. SolverInitializeK3biilu2 (" << file_options << ") \n";)//db!
    if (file_options == NULL) return;
    std::string s = file_options;
    if (s == "" || s == " ") return;
    std::ifstream is;
    T(std::cout << "##### ins. SolverInitializeK3biilu2: bef. open(" << file_options << ") \n";)//db!
    is.open(file_options, std::ifstream::in);
    k3d::SParams parPrec = *(data->pParams);
    ParIter parIter = *(data->pParIter);
    if (s == "k3biilu2_options.txt") {
        getline(is, s);
        sscanf(s.c_str(), "%d", &parPrec.prec_float);   //01 prec_float
        getline(is, s);
        sscanf(s.c_str(), "%d", &parPrec.ncycle);       //02 ncycle
        getline(is, s);
        sscanf(s.c_str(), "%d", &parPrec.ordtype);      //03 ordtype
        getline(is, s);
        sscanf(s.c_str(), "%d", &parPrec.collap);       //04 collap
        getline(is, s);
        sscanf(s.c_str(), "%d", &parPrec.sctype);       //05 sctype
        getline(is, s);
        sscanf(s.c_str(), "%d", &parPrec.nitersc);      //06 nitersc
        getline(is, s);
        sscanf(s.c_str(), "%d", &parPrec.fcttype);      //07 fcttype
        getline(is, s);
        sscanf(s.c_str(), "%lg", &parPrec.pivmin);       //08 pivmin
        getline(is, s);
        sscanf(s.c_str(), "%lg", &parPrec.tau1);         //09 tau1
        getline(is, s);
        sscanf(s.c_str(), "%lg", &parPrec.tau2);         //10 tau2
        getline(is, s);
        sscanf(s.c_str(), "%lg", &parPrec.theta);        //11 theta
        getline(is, s);
        sscanf(s.c_str(), "%d", &parIter.ittype);       //12 ittype
        getline(is, s);
        sscanf(s.c_str(), "%d", &parIter.niter_cycle);  //13 niter_cycle
        getline(is, s);
        sscanf(s.c_str(), "%d", &parIter.ncoef);        //14 ncoef
        getline(is, s);
        sscanf(s.c_str(), "%d", &parIter.niter_cycle2); //15 niter_cycle2
        getline(is, s);
        sscanf(s.c_str(), "%lg", &parIter.eps);          //16 eps
        getline(is, s);
        sscanf(s.c_str(), "%d", &parIter.maxit);        //17 maxit
        getline(is, s);
        sscanf(s.c_str(), "%d", &parIter.ichk);         //18 ichk
        getline(is, s);
        sscanf(s.c_str(), "%d", &parIter.msglev);       //19 msglev
        data->parameters_initialized = true;
    } else if (s == "ctrl_dat") {
        getline(is, s);                                                  //1 skip iext
        getline(is, s);                                                  //2 skip mtx filename
        getline(is, s);                                                  //3 skip rhs filename
        getline(is, s);
        sscanf(s.c_str(), "%lg", &parIter.eps);          //4 eps
        getline(is, s);
        sscanf(s.c_str(), "%d", &parIter.maxit);        //5 maxit
        getline(is, s);
        sscanf(s.c_str(), "%d", &parIter.niter_cycle);  //6 kgmr
        getline(is, s);
        sscanf(s.c_str(), "%d", &parIter.ncoef);        //7 kdeg
        getline(is, s);
        sscanf(s.c_str(), "%d", &parPrec.ncycle);       //8 kovl
        getline(is, s);
        sscanf(s.c_str(), "%lg", &parPrec.tau1);         //9 tau
        parPrec.tau2 = -1.0;
        parIter.ittype = (parIter.niter_cycle > 1) ? 1 : 0;
        data->parameters_initialized = true;
        //? msglev
    } else if (s == "biilu2_options.txt") { // file: "biilu2_options.txt"
        getline(is, s);
        sscanf(s.c_str(), "%d", &parPrec.ncycle); //1 kovl
        getline(is, s);
        sscanf(s.c_str(), "%lg", &parPrec.tau1);   //2 tau
        getline(is, s);
        sscanf(s.c_str(), "%lg", &parIter.eps);    //3 eps
        getline(is, s);
        sscanf(s.c_str(), "%d", &parIter.maxit);  //4 maxit
        getline(is, s);
        sscanf(s.c_str(), "%d", &parIter.msglev); //5 msglev
        parPrec.tau2 = -1.0;
        data->parameters_initialized = true;
    }
    T(std::cout << "##### ins. SolverInitializeK3biilu2:  prec_float=" << parPrec.prec_float << " ncycle=" << parPrec.ncycle << " ordtype=" << parPrec.ordtype
                << " collap=" << parPrec.collap << " sctype=" << parPrec.sctype << " nitersc=" << parPrec.nitersc << " fcttype=" << parPrec.fcttype
                << " pivmin=" << parPrec.pivmin << " tau1=" << parPrec.tau1 << " tau2=" << parPrec.tau2 << " theta=" << parPrec.theta << " ittype="
                << parIter.ittype << " niter_cycle=" << parIter.niter_cycle << " ncoef=" << parIter.ncoef << " niter_cycle2=" << parIter.niter_cycle2 << " eps="
                << parIter.eps << " maxit=" << parIter.maxit << " ichk=" << parIter.ichk << " msglev=" << parIter.msglev << " from: " << file_options
                << " \n";)//db!
    (void) argc;
    (void) argv;
}

bool SolverIsFinalizedK3biilu2() {
    return true; //no need to finalize
}

void SolverFinalizeK3biilu2() {
}

void SolverDestroyDataK3biilu2(bcg_k3biilu2 **data) {
    if (data != NULL) {
        if (*data != NULL) {
            freebcg_k3(*data);
            free(*data);
        }
        *data = NULL;
    }
}

void SolverInitDataK3biilu2(bcg_k3biilu2 **data, INMOST_MPI_Comm comm, const char *name) {
    T(std::cout << "##### ins. SolverInitDataK3biilu2 \n";)//db!
    *data = (bcg_k3biilu2 *) malloc(sizeof(bcg_k3biilu2));
    (*data)->n = 0;
    (*data)->nproc = 0;
    (*data)->pSolver = new k3d::CK3D_Solver<int, double, double>();
    (*data)->pParIter = new ParIter();
    (*data)->pParams = new k3d::SParams();
    (*data)->parameters_initialized = false;
    ParametersDefault(*((*data)->pParIter), *((*data)->pParams));
    T(cout << "HHHHH bef. Clean in SolverInitDataK3biilu2 \n";)//DB!
    (*data)->pSolver->Clean();
    (void) comm;
    (void) name;
}

void SolverCopyDataK3biilu2(bcg_k3biilu2 **data, bcg_k3biilu2 *other_data, INMOST_MPI_Comm comm) {
    throw INMOST::NotImplemented; //later
    (void) data;
    (void) other_data;
    (void) comm;
}

void SolverAssignDataK3biilu2(bcg_k3biilu2 *data, bcg_k3biilu2 *other_data) {
    throw INMOST::NotImplemented; //later
    (void) data;
    (void) other_data;
}

void SolverSetMatrixK3biilu2(bcg_k3biilu2 *data, matrix_k3biilu2 *matrix_data, bool same_pattern, bool reuse_preconditioner) {
    T(std::cout << "##### ins. SolverSetMatrixK3biilu2 \n";)//db!
    //if( A == NULL) std::cout<<"##### A == NULL ... \n";//db!
    //if( m == NULL) std::cout<<"##### m == NULL ... \n";//db!
    if (matrix_data == NULL || data == NULL) throw INMOST::DataCorruptedInSolver;
    T(std::cout << "##### ins. SolverSetMatrixK3biilu2 bef. initbcg or newmatrixbcg \n";)//db!
    if (data->n == 0)
        initbcg_k3(data, matrix_data, data->pParIter->eps);
    else
        newmatrixbcg_k3(data, matrix_data, reuse_preconditioner);
    (void) same_pattern;
    T(std::cout << "##### ins. SolverSetMatrixK3biilu2 bef. return \n";)//db!
}

bool SolverSolveK3biilu2(bcg_k3biilu2 *data, vector_k3biilu2 *rhs_data, vector_k3biilu2 *sol_data) {
    T(std::cout << "##### ins. SolverSolveK3biilu2 \n";)//db!
    return solvebcg_k3(data, rhs_data, sol_data) == 0;
}
