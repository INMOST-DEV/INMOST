#include "solver_k3biilu2.h"

//#if defined(USE_SOLVER_K3BIILU2)

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <string>

#include "k3d.h"


static int    set_kovl   = 0;    // number of overlap layers: kovl=0,1,2,...
static double set_tau    = 3e-3; // the ILU2 precision (for the submatrix factorization); tau=3e-3
static double set_eps    = 1e-5; // the residual precision: ||r|| < eps * ||b||; eps=1e-6
static int    set_nit    = 999;  // number of iterations permitted; nit=999
static int    set_msglev = 2;    // messages level; msglev=0 for silent; msglev=1 to output solution statistics

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
} Vector;

/*****************************************************************************/
#include <stdlib.h>   // for malloc()
#include <stdio.h>    // for printf()
#if defined (USE_MPI)
#include <mpi.h>      // for MPI_COMM_WORLD etc.
#endif

#if defined (USE_MPI)
#define MPI_BARRIER(comm) \
        MPI_Barrier(comm)
#else
#define MPI_BARRIER(comm)
#endif

#if defined (USE_MPI)
#define MPI_ALLREDUCE(sendbuf, recvbuf, count, datatype, op, comm) \
        MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm)
#else
#define MPI_ALLREDUCE(sendbuf, recvbuf, count, datatype, op, comm)
#endif

// Solve the linear system A X = B by
//  k3biilu2_bcg solver with working memory allocations and statistics output
int k3biilu2_bcg (
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
        int    *istat,  //[16],  // integer statistics array on return
        double *dstat)  //[16]); // double  statistics array on return
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
{
    // Initialize MPI variables
    int np=1, mp=0;
    MPI_Comm comm = MPI_COMM_WORLD; //TODO: MSPP_COMM_WORLD
#if defined (USE_MPI)
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &mp);
#endif

    k3d::SSolverParams params;

    params.prec_float = 1;
    params.ncycle = kovl;
    params.fcttype = -1;
    params.tau1 = tau;
    params.tau2 = tau*tau;
    params.theta = 0.1e0;

    std::vector<long long> blks (np+1);
    std::vector<int> blk2cpu (np+1);

    long long *pblks = &blks[0];
    int *pblk2cpu = &blk2cpu[0];

    int i;

    for (i=0;i<=np;i++) pblks[i] = ibl[i];
    for (i=0;i<np;i++) pblk2cpu[i] = i;

    k3d::CSolver<int,double,double> *pSolver = new k3d::CSolver<int,double,double>;

    pSolver->PrepareSolver ((void *)&comm, &params, np, pblks, pblk2cpu,
                            true, ia, ja, a);

    if (nit > 0) {
       int ichk = 5;
       msglev = 0;
       ofstream *pfout = NULL;
       if (mp == 0) msglev = 2;
       double rhs_norm, res_ini, res_fin;
       int niter;
       pSolver->BiCGStab (nit, eps, ichk, msglev, pfout,
                          b, x,
			  rhs_norm, res_ini, niter, res_fin);
       istat[2] = niter;
       dstat[2] = (rhs_norm == 0e0) ? res_fin : res_fin/rhs_norm;
    }
    
    *ierr = 0;

    return *ierr;
}
/*****************************************************************************/

/* Initialize bcg solver */
static int initbcg(bcg *s, matrix *A, double eps);
/* Reinitialize solver preconditioner with new matrix A */
static int renewbcg(bcg *s, double *A);
/* Solve linear system */
/*static*/ int solvebcg(bcg *s, Vector *b, Vector *x);
/* Free memory used by solver */
static void freebcg(bcg *s);

/*****************************************************************************/

/* Initialize solver with new matrix A */
static int newmatrixbcg(bcg *s, matrix *A, bool same_precond) 
{
    if( s->n != A->n && same_precond ) throw INMOST::CannotReusePreconditionerOfDifferentSize;
    s->n     = A->n;
    s->nproc = A->nproc;
    s->ibl   = A->ibl;
    s->ia    = A->ia;
    s->ja    = A->ja;
    if( !same_precond )
    {
	//do nothing...
        //std::cout<<"##### inside newmatrixbcg bef. renewbcg \n";//db!
	return renewbcg(s, A->A);
    }
    else return 0;
}

/* solver */
/* Initialize bcg solver */
int initbcg(bcg *s, matrix *A, double eps) 
{
    s->eps = eps;
    //std::cout<<"##### inside initbcg bef. newmatrixbcg eps="<<eps<<" \n";//db!
    return newmatrixbcg(s, A, false);
}

/* Reinitialize solver preconditioner with new matrix A */
int renewbcg(bcg *s, double *A) 
{
    //reinitialize matrix values
    s->a = A;
    //BIILU2 preconditioner construction...
    s->kovl    = set_kovl;
    s->tau     = set_tau;
    s->msglev  = set_msglev;

    int   job  = 0;
    int maxit  = 0;
    s->len_r8  = 0; //to be set by nbl inside k3biilu2_bcg
    if(s->W) free(s->W); s->W=NULL;

    int ierr = 0;
    //std::cout<<"##### inside renewbcg bef. k3biilu2_bcg\n";//db!
    //double *B = (double*) malloc(sizeof(double)*s->n); //db!!!!!!!!!!!!!!
    //double *X = (double*) malloc(sizeof(double)*s->n); //db!!!!!!!!!!!!!!
    k3biilu2_bcg (s->ibl, s->ia, s->ja, s->a,
                NULL, NULL,
		job, &s->len_r8, &s->W,
                s->kovl, s->tau, s->eps, maxit, s->msglev,
	        &ierr, s->istat, s->dstat);
    //std::cout<<"##### inside renewbcg aft. k3biilu2_bcg\n";//db!
    //free(B);free(X);//db!!!!!!!!!!!!!!!
    if (ierr) printf("initialization of biilu2 failed, ierr=%d\n", ierr);

    return ierr;
}

/* Solve linear system */
int solvebcg(bcg *s, Vector *b, Vector *x) 
{
    s->kovl    = set_kovl;
    s->tau     = set_tau;
//  s->eps     = set_eps;
    s->nit     = set_nit;
    s->msglev  = set_msglev;

    int   job  = 1;
    int maxit  = s->nit;

    int ierr = 0;
    //std::cout<<"##### inside solvebcg bef. k3biilu2_bcg\n";//db!
    k3biilu2_bcg (s->ibl, s->ia, s->ja, s->a, b->v, x->v,
		job, &s->len_r8, &s->W,
                s->kovl, s->tau, s->eps, maxit, s->msglev,
	        &ierr, s->istat, s->dstat);
    //std::cout<<"##### inside solvebcg aft. k3biilu2_bcg\n";//db!

    s->ITER  = s->istat[2];
    s->RESID = s->dstat[2];

    return ierr;
}

/* Free memory used by solver */
void freebcg(bcg *s) 
{
    if(s->W) free(s->W); s->W=NULL; //-IK!!!!!!!!!!!!!!!!!!
}

/*****************************************************************************/

void MatrixCopyDataK3biilu2(void ** ppA, void * pB)
{
	matrix * B = (matrix *)pB;
	if( ppA == NULL || pB == NULL ) throw INMOST::DataCorruptedInSolver;
	*ppA = malloc(sizeof(matrix));
	matrix * A = (matrix *)ppA;
	A->n = B->n;
	if( B->n != 0 )
	{
		int nnz = B->ia[B->n] - B->ia[0];
		A->nproc = B->nproc;
		A->ibl = (int *) malloc(sizeof(int)*(A->nproc+1));
		memcpy(A->ibl,B->ibl,sizeof(int)*(A->nproc+1));
		A->ia = (int *) malloc(sizeof(int)*(A->n+1));
		memcpy(A->ia,B->ia,sizeof(int)*(A->n+1));
		A->ja = (int *) malloc(sizeof(int)*nnz);
		memcpy(A->ja,B->ja,sizeof(int)*nnz);
		A->A = (double *) malloc(sizeof(double)*nnz);
		memcpy(A->A,B->A,sizeof(double)*nnz);
	}
}

void MatrixAssignDataK3biilu2(void * pA, void * pB)
{
	matrix * A = (matrix *)pA;
	matrix * B = (matrix *)pB;
	if( A == NULL || B == NULL ) throw INMOST::DataCorruptedInSolver;
	if( A != B )
	{
		if( A->n != 0 )
		{
			free(A->ibl);
			free(A->ia);
			free(A->ja);
			free(A->A);
		}
		if( B->n != 0 )
		{
			int nnz = B->ia[B->n] - B->ia[0];
			A->n = B->n;
			A->nproc = B->nproc;
			A->ibl = (int *) malloc(sizeof(int)*(A->nproc+1));
			memcpy(A->ibl,B->ibl,sizeof(int)*(A->nproc+1));
			A->ia = (int *) malloc(sizeof(int)*(A->n+1));
			memcpy(A->ia,B->ia,sizeof(int)*(A->n+1));
			A->ja = (int *) malloc(sizeof(int)*nnz);
			memcpy(A->ja,B->ja,sizeof(int)*nnz);
			A->A = (double *) malloc(sizeof(double)*nnz);
			memcpy(A->A,B->A,sizeof(double)*nnz);	
		}
	}
}

void MatrixInitDataK3biilu2(void ** ppA, INMOST_MPI_Comm comm, const char * name)
{
        //std::cout<<"##### ins. MatrixInitDataK3biilu2 \n";//db!
	if( ppA == NULL ) throw INMOST::DataCorruptedInSolver;
	if( *ppA == NULL )
	{
		*ppA = malloc(sizeof(matrix));
		matrix * A = (matrix *)*ppA;
		A->n = 0;
		A->nproc = 0;
                //std::cout<<"##### ins. MatrixInitDataK3biilu2 n=nproc=0 \n";//db!
	}
    (void) comm;
    (void) name;
}

void MatrixDestroyDataK3biilu2(void ** pA)
{
	matrix * A = (matrix *)(*pA);
	if( A != NULL )
	{
		if( A->n != 0 )
		{
			free(A->ibl);
			free(A->ia);
			free(A->ja);
			free(A->A);
                        //std::cout<<"##### ins. MatrixDestroyDataK3biilu2 ...free \n";//db!
		}
		free(*pA);
		*pA = NULL;
	}
}



void MatrixFillK3biilu2(void * pA, int size, int nproc, int * ibl, int * ia, int * ja, double * values)
{
        //std::cout<<"##### ins. MatrixFillK3biilu2 n="<<size<<" nproc="<<nproc<<" \n";//db!
	if( pA == NULL ) throw INMOST::DataCorruptedInSolver;
	matrix * A = (matrix *) pA;
	A->n = size;
	A->nproc = nproc;
	A->ibl = ibl;
	A->ia = ia;
	A->ja = ja;
	A->A = values;
}

void MatrixFillValuesK3biilu2(void * pA, double * values)
{
        //std::cout<<"##### ins. MatrixFillValuesK3biilu2 \n";//db!
	if( pA == NULL ) throw INMOST::DataCorruptedInSolver;
	matrix * A = (matrix *) pA;
	free(A->A);
	A->A = values;
}

void MatrixFinalizeK3biilu2(void * data)
{
	//don't need to do anything
    (void) data;
}

void VectorInitDataK3biilu2(void ** ppA, INMOST_MPI_Comm comm, const char * name)
{
	if( ppA == NULL ) throw INMOST::DataCorruptedInSolver;
	*ppA = malloc(sizeof(Vector));
	Vector * A = (Vector *)*ppA;
	A->n = 0;
    (void) comm;
    (void) name;
}

void VectorCopyDataK3biilu2(void ** ppA, void * pB)
{
        //std::cout<<"##### ins. VectorCopyDataK3biilu2 \n";//db!
	if( ppA == NULL || pB == NULL ) throw INMOST::DataCorruptedInSolver;
	*ppA = malloc(sizeof(Vector));
	Vector * A = (Vector *)*ppA;
	Vector * B = (Vector *)pB;
	A->n = B->n;
	if( B->n != 0 )
	{
		A->v = (double *)malloc(sizeof(double)*A->n);
		memcpy(A->v,B->v,sizeof(double)*A->n);
	}
}

void VectorAssignDataK3biilu2(void * pA, void * pB)
{
        //std::cout<<"##### ins. VectorAssignDataK3biilu2 \n";//db!
	Vector * A = (Vector *)pA;
	Vector * B = (Vector *)pB;
	if( A == NULL || B == NULL ) throw INMOST::DataCorruptedInSolver;
	if( A != B )
	{
		if( A->n != 0 ) free(A->v);
		A->n = B->n;
		if( B->n != 0 )
		{
			A->v = (double *) malloc(sizeof(double)*A->n);
			memcpy(A->v,B->v,sizeof(double)*A->n);
		}
	}
}

void VectorPreallocateK3biilu2(void * pA, int size)
{
	Vector * A = (Vector *)pA;
	if( A == NULL ) throw INMOST::DataCorruptedInSolver;
	A->n = size;
	A->v = (double *)malloc(sizeof(double)*size);
}

void VectorFillK3biilu2(void * pA, double * values)
{
	Vector * A = (Vector *)pA;
	if( A == NULL ) throw INMOST::DataCorruptedInSolver;
	memcpy(A->v,values,sizeof(double)*A->n);
}
void VectorLoadK3biilu2(void * pA, double * values)
{
	Vector * A = (Vector *)pA;
	if( A == NULL ) throw INMOST::DataCorruptedInSolver;
	memcpy(values,A->v,sizeof(double)*A->n);
}

void VectorFinalizeK3biilu2(void * data)
{
    (void) data;
}

void VectorDestroyDataK3biilu2(void ** ppA)
{
	if( ppA == NULL) throw INMOST::DataCorruptedInSolver;
	if( *ppA != NULL )
	{
		Vector * A = (Vector *)*ppA;
		free(A->v);
		free(*ppA);
		*ppA = NULL;
	}
}

void SolverInitializeK3biilu2(int * argc, char *** argv, const char * file_options)
{
    //std::cout<<"##### ins. SolverInitializeK3biilu2 ("<<file_options<<") \n";//db!
    if (file_options == NULL) return;
    std::string s = file_options;
    if (s == "" || s == " ") return;
    std::ifstream is;
    //std::cout<<"##### ins. SolverInitializeK3biilu2: bef. open("<<file_options<<") \n";//db!
    is.open(file_options, std::ifstream::in);
    if (s == "ctrl_dat") {
        getline(is, s);                                      //1 skip iext
        getline(is, s);                                      //2 skip mtx filename
        getline(is, s);                                      //3 skip rhs filename
        getline(is, s); sscanf(s.c_str(), "%lg", &set_eps);  //4 eps
        getline(is, s); sscanf(s.c_str(), "%d",  &set_nit);  //5 nit
        getline(is, s); sscanf(s.c_str(), "%d",  &set_kovl); //6 kovl
        getline(is, s); sscanf(s.c_str(), "%lg", &set_tau);  //7 tau
        //? msglev
    } else {
        getline(is, s); sscanf(s.c_str(), "%d",  &set_kovl);   //2 kovl
        getline(is, s); sscanf(s.c_str(), "%lg", &set_tau);    //3 tau
        getline(is, s); sscanf(s.c_str(), "%lg", &set_eps);    //4 eps
        getline(is, s); sscanf(s.c_str(), "%d",  &set_nit);    //5 nit
        getline(is, s); sscanf(s.c_str(), "%d",  &set_msglev); //6 msglev
    }
    //std::cout<<"##### ins. SolverInitializeK3biilu2:  kovl="<<set_kovl<<" tau="<<set_tau<<" eps="<<set_eps<<" nit="<<set_nit<<" msglev="<<set_msglev<<" \n";//db!
	(void) argc;
	(void) argv;
}

bool SolverIsFinalizedK3biilu2()
{
	return true; //no need to finalize
}

void SolverFinalizeK3biilu2()
{
}

void SolverDestroyDataK3biilu2(void ** data)
{
	if( data != NULL )
	{
		if( *data != NULL )
		{
			bcg * m = (bcg *)*data;
			freebcg(m);
			free(m);
		}
		*data = NULL;
	}
}

void SolverInitDataK3biilu2(void ** data, INMOST_MPI_Comm comm, const char * name)
{
        //std::cout<<"##### ins. SolverInitDataK3biilu2 \n";//db!
	*data = malloc(sizeof(bcg));
	((bcg *)*data)->n = 0;
	((bcg *)*data)->nproc = 0;
	((bcg *)*data)->len_r8 = 0;
	((bcg *)*data)->W = NULL;
	(void) comm;
	(void) name;
}

void SolverCopyDataK3biilu2(void ** data, void * other_data, INMOST_MPI_Comm comm)
{
	throw INMOST::NotImplemented; //later
	(void) data;
	(void) other_data;
	(void) comm;
}

void SolverAssignDataK3biilu2(void * data, void * other_data)
{
	throw INMOST::NotImplemented; //later
	(void) data;
	(void) other_data;
}

void SolverSetMatrixK3biilu2(void * data, void * matrix_data, bool same_pattern, bool reuse_preconditioner)
{
        //std::cout<<"##### ins. SolverSetMatrixK3biilu2 \n";//db!
	bcg * m = (bcg *)data;
	matrix * A = (matrix *)matrix_data;
        //if( A == NULL) std::cout<<"##### A == NULL ... \n";//db!
        //if( m == NULL) std::cout<<"##### m == NULL ... \n";//db!
	if( A == NULL || m == NULL ) throw INMOST::DataCorruptedInSolver;
        //std::cout<<"##### ins. SolverSetMatrixK3biilu2 bef. initbcg or newmatrixbcg \n";//db!
	if( m->n == 0 )
		initbcg(m,A,set_eps);
	else
		newmatrixbcg(m,A,reuse_preconditioner);
	(void) same_pattern;
        //std::cout<<"##### ins. SolverSetMatrixK3biilu2 bef. return \n";//db!
}

bool SolverSolveK3biilu2(void * data, void * rhs_data, void * sol_data)
{
        //std::cout<<"##### ins. SolverSolveK3biilu2 \n";//db!
	bcg * m = (bcg *)data;
	Vector * rhs = (Vector*)rhs_data, * sol = (Vector *)sol_data;
	return solvebcg(m,rhs,sol) == 0;
}

int SolverIterationNumberK3biilu2(void * data)
{
	return ((bcg *)data)->ITER;
}

double SolverResidualNormK3biilu2(void * data)
{
	return ((bcg *)data)->RESID;
}

/*
void SolverAddOtherStatK3biilu2(void * data, unsigned int * pivmod, double * prdens, double * t_prec, double * t_iter)
{
	*pivmod += ((bcg *)data)->istat[0];
	*prdens += ((bcg *)data)->dstat[0];
	*t_prec += ((bcg *)data)->dstat[7];
	*t_iter += ((bcg *)data)->dstat[9];
	return;
}
*/

//#endif //USE_SOLVER_K3BIILU2
