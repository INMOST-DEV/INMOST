//#include "inmost_solver.h"
//#if defined(USE_SOLVER)
//#if defined(USE_SOLVER_FCBIILU2)

#include "solver_fcbiilu2.h"

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <string>


static int    set_kovl   = 0;    // number of overlap layers: kovl=0,1,2,...
static double set_tau    = 3e-3; // the ILU2 precision (for the submatrix factorization); tau=3e-3
static double set_eps    = 1e-5; // the residual precision: ||r|| < eps * ||b||; eps=1e-6
static int    set_nit    = 999;  // number of iterations permitted; nit=999
//static int    set_kgmr   = 1;    // iterative scheme (BiCGStab: kgmr=1;  GMR[kdeg]: kgmr>3*kdeg+2)
//static int    set_kdeg   = 1;    // degree of polynomial for GMR[d] (ignored if kgmr=1)
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
} vector;

/*****************************************************************************/

//----- File: biilu2.h if available
// Solve the linear system A X = B by
//  biilu2_bcg solver with working memory allocations and statistics output
int biilu2_bcg (
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
static int initbcg(bcg *s, matrix *A, double eps);
/* Reinitialize solver preconditioner with new matrix A */
static int renewbcg(bcg *s, double *A);
/* Solve linear system */
/*static*/ int solvebcg(bcg *s, vector *b, vector *x);
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
    s->len_r8  = 0; //to be set by nbl inside biilu2_bcg
    if(s->W) free(s->W); s->W=NULL;

    int ierr = 0;
    //std::cout<<"##### inside renewbcg bef. biilu2_bcg\n";//db!
    double *B = (double*) malloc(sizeof(double)*s->n); //db!!!!!!!!!!!!!!
    //double *X = (double*) malloc(sizeof(double)*s->n); //db!!!!!!!!!!!!!!
    biilu2_bcg (s->ibl, s->ia, s->ja, s->a,
                //B, X, //!!!!!KAPORIN!!!!!!!!!!!!!
                B, NULL, //!!!!!KAPORIN!!!!!!!!!!!!
                //NULL, NULL, //!!!!!INK!!!!!!!!!!!
		job, &s->len_r8, &s->W,
                s->kovl, s->tau, s->eps, maxit, s->msglev,
	        &ierr, s->istat, s->dstat);
    //std::cout<<"##### inside renewbcg aft. biilu2_bcg\n";//db!
    free(B); //free(X);//db!!!!!!!!!!!!!!!
    if (ierr) printf("initialization of biilu2 failed, ierr=%d\n", ierr);

    return ierr;
}

/* Solve linear system */
int solvebcg(bcg *s, vector *b, vector *x) 
{
    s->kovl    = set_kovl;
    s->tau     = set_tau;
//  s->eps     = set_eps;
    s->nit     = set_nit;
    s->msglev  = set_msglev;

    int   job  = 1;
    int maxit  = s->nit;

    int ierr = 0;
    //std::cout<<"##### inside solvebcg bef. biilu2_bcg\n";//db!
    biilu2_bcg (s->ibl, s->ia, s->ja, s->a, b->v, x->v,
		job, &s->len_r8, &s->W,
                s->kovl, s->tau, s->eps, maxit, s->msglev,
	        &ierr, s->istat, s->dstat);
    //std::cout<<"##### inside solvebcg aft. biilu2_bcg\n";//db!

    s->ITER  = s->istat[2];
    s->RESID = s->dstat[2];

    return ierr;
}

/* Free memory used by solver */
void freebcg(bcg *s) 
{
    if(s->W) free(s->W); s->W=NULL;
}

/*****************************************************************************/

void MatrixCopyDataFcbiilu2(void ** ppA, void * pB)
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

void MatrixAssignDataFcbiilu2(void * pA, void * pB)
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

void MatrixInitDataFcbiilu2(void ** ppA, INMOST_MPI_Comm comm, const char * name)
{
        //std::cout<<"##### ins. MatrixInitDataFcbiilu2 \n";//db!
	if( ppA == NULL ) throw INMOST::DataCorruptedInSolver;
	if( *ppA == NULL )
	{
		*ppA = malloc(sizeof(matrix));
		matrix * A = (matrix *)*ppA;
		A->n = 0;
		A->nproc = 0;
                //std::cout<<"##### ins. MatrixInitDataFcbiilu2 n=nproc=0 \n";//db!
	}
    (void) comm;
    (void) name;
}

void MatrixDestroyDataFcbiilu2(void ** pA)
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
                        //std::cout<<"##### ins. MatrixDestroyDataFcbiilu2 ...free \n";//db!
		}
		free(*pA);
		*pA = NULL;
	}
}



void MatrixFillFcbiilu2(void * pA, int size, int nproc, int * ibl, int * ia, int * ja, double * values)
{
        //std::cout<<"##### ins. MatrixFillFcbiilu2 n="<<size<<" nproc="<<nproc<<" \n";//db!
	if( pA == NULL ) throw INMOST::DataCorruptedInSolver;
	matrix * A = (matrix *) pA;
	A->n = size;
	A->nproc = nproc;
	A->ibl = ibl;
	A->ia = ia;
	A->ja = ja;
	A->A = values;
}

void MatrixFillValuesFcbiilu2(void * pA, double * values)
{
        //std::cout<<"##### ins. MatrixFillValuesFcbiilu2 \n";//db!
	if( pA == NULL ) throw INMOST::DataCorruptedInSolver;
	matrix * A = (matrix *) pA;
	free(A->A);
	A->A = values;
}

void MatrixFinalizeFcbiilu2(void * data)
{
	//don't need to do anything
    (void) data;
}

void VectorInitDataFcbiilu2(void ** ppA, INMOST_MPI_Comm comm, const char * name)
{
	if( ppA == NULL ) throw INMOST::DataCorruptedInSolver;
	*ppA = malloc(sizeof(vector));
	vector * A = (vector *)*ppA;
	A->n = 0;
    (void) comm;
    (void) name;
}

void VectorCopyDataFcbiilu2(void ** ppA, void * pB)
{
        //std::cout<<"##### ins. VectorCopyDataFcbiilu2 \n";//db!
	if( ppA == NULL || pB == NULL ) throw INMOST::DataCorruptedInSolver;
	*ppA = malloc(sizeof(vector));
	vector * A = (vector *)*ppA;
	vector * B = (vector *)pB;
	A->n = B->n;
	if( B->n != 0 )
	{
		A->v = (double *)malloc(sizeof(double)*A->n);
		memcpy(A->v,B->v,sizeof(double)*A->n);
	}
}

void VectorAssignDataFcbiilu2(void * pA, void * pB)
{
        //std::cout<<"##### ins. VectorAssignDataFcbiilu2 \n";//db!
	vector * A = (vector *)pA;
	vector * B = (vector *)pB;
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

void VectorPreallocateFcbiilu2(void * pA, int size)
{
	vector * A = (vector *)pA;
	if( A == NULL ) throw INMOST::DataCorruptedInSolver;
	A->n = size;
	A->v = (double *)malloc(sizeof(double)*size);
}

void VectorFillFcbiilu2(void * pA, double * values)
{
	vector * A = (vector *)pA;
	if( A == NULL ) throw INMOST::DataCorruptedInSolver;
	memcpy(A->v,values,sizeof(double)*A->n);
}
void VectorLoadFcbiilu2(void * pA, double * values)
{
	vector * A = (vector *)pA;
	if( A == NULL ) throw INMOST::DataCorruptedInSolver;
	memcpy(values,A->v,sizeof(double)*A->n);
}

void VectorFinalizeFcbiilu2(void * data)
{
    (void) data;
}

void VectorDestroyDataFcbiilu2(void ** ppA)
{
	if( ppA == NULL) throw INMOST::DataCorruptedInSolver;
	if( *ppA != NULL )
	{
		vector * A = (vector *)*ppA;
		free(A->v);
		free(*ppA);
		*ppA = NULL;
	}
}

void SolverInitializeFcbiilu2(int * argc, char *** argv, const char * file_options)
{
    //std::cout<<"##### ins. SolverInitializeFcbiilu2 ("<<file_options<<") \n";//db!
    if (file_options == NULL) return;
    std::string s = file_options;
    if (s == "" || s == " ") return;
    std::ifstream is;
    //std::cout<<"##### ins. SolverInitializeFcbiilu2: bef. open("<<file_options<<") \n";//db!
    is.open(file_options, std::ifstream::in);
    if (s == "ctrl_dat") {
        getline(is, s);                                      //1 skip ipart
        getline(is, s);                                      //2 skip mtx filename
        getline(is, s);                                      //3 skip rhs filename
        getline(is, s); sscanf(s.c_str(), "%lg", &set_eps);  //4 eps
        getline(is, s); sscanf(s.c_str(), "%d",  &set_nit);  //5 nit
        getline(is, s); //sscanf(s.c_str(), "%d",  &set_kgmr); //6- skip kgmr
        getline(is, s); //sscanf(s.c_str(), "%d",  &set_kdeg); //7- skip kdeg
        getline(is, s); sscanf(s.c_str(), "%d",  &set_kovl); //8 kovl
        getline(is, s); sscanf(s.c_str(), "%lg", &set_tau);  //9 tau
        //? msglev
    } else { // file: "biilu2_options.txt"
        getline(is, s); sscanf(s.c_str(), "%d",  &set_kovl);   //1 kovl
        getline(is, s); sscanf(s.c_str(), "%lg", &set_tau);    //2 tau
        getline(is, s); sscanf(s.c_str(), "%lg", &set_eps);    //3 eps
        getline(is, s); sscanf(s.c_str(), "%d",  &set_nit);    //4 nit
        getline(is, s); sscanf(s.c_str(), "%d",  &set_msglev); //5 msglev
    }
    //std::cout<<"##### ins. SolverInitializeFcbiilu2:  kovl="<<set_kovl<<" tau="<<set_tau<<" eps="<<set_eps<<" nit="<<set_nit<<" msglev="<<set_msglev<<" from: "<<file_options<<" \n";//db!
	(void) argc;
	(void) argv;
}

bool SolverIsFinalizedFcbiilu2()
{
	return true; //no need to finalize
}

void SolverFinalizeFcbiilu2()
{
}

void SolverDestroyDataFcbiilu2(void ** data)
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

void SolverInitDataFcbiilu2(void ** data, INMOST_MPI_Comm comm, const char * name)
{
        //std::cout<<"##### ins. SolverInitDataFcbiilu2 \n";//db!
	*data = malloc(sizeof(bcg));
	((bcg *)*data)->n = 0;
	((bcg *)*data)->nproc = 0;
	((bcg *)*data)->len_r8 = 0;
	((bcg *)*data)->W = NULL;
	(void) comm;
	(void) name;
}

void SolverCopyDataFcbiilu2(void ** data, void * other_data, INMOST_MPI_Comm comm)
{
	throw INMOST::NotImplemented; //later
	(void) data;
	(void) other_data;
	(void) comm;
}

void SolverAssignDataFcbiilu2(void * data, void * other_data)
{
	throw INMOST::NotImplemented; //later
	(void) data;
	(void) other_data;
}

void SolverSetMatrixFcbiilu2(void * data, void * matrix_data, bool same_pattern, bool reuse_preconditioner)
{
        //std::cout<<"##### ins. SolverSetMatrixFcbiilu2 \n";//db!
	bcg * m = (bcg *)data;
	matrix * A = (matrix *)matrix_data;
        //if( A == NULL) std::cout<<"##### A == NULL ... \n";//db!
        //if( m == NULL) std::cout<<"##### m == NULL ... \n";//db!
	if( A == NULL || m == NULL ) throw INMOST::DataCorruptedInSolver;
        //std::cout<<"##### ins. SolverSetMatrixFcbiilu2 bef. initbcg or newmatrixbcg \n";//db!
	if( m->n == 0 )
		initbcg(m,A,set_eps);
	else
		newmatrixbcg(m,A,reuse_preconditioner);
	(void) same_pattern;
        //std::cout<<"##### ins. SolverSetMatrixFcbiilu2 bef. return \n";//db!
}

bool SolverSolveFcbiilu2(void * data, void * rhs_data, void * sol_data)
{
        //std::cout<<"##### ins. SolverSolveFcbiilu2 \n";//db!
	bcg * m = (bcg *)data;
	vector * rhs = (vector*)rhs_data, * sol = (vector *)sol_data;
	return solvebcg(m,rhs,sol) == 0;
}

int SolverIterationNumberFcbiilu2(void * data)
{
	return ((bcg *)data)->ITER;
}

double SolverResidualNormFcbiilu2(void * data)
{
	return ((bcg *)data)->RESID;
}

/*
void SolverAddOtherStatFcbiilu2(void * data, unsigned int * pivmod, double * prdens, double * t_prec, double * t_iter)
{
	*pivmod += ((bcg *)data)->istat[0];
	*prdens += ((bcg *)data)->dstat[0];
	*t_prec += ((bcg *)data)->dstat[7];
	*t_iter += ((bcg *)data)->dstat[9];
	return;
}
*/

//#endif //USE_SOLVER
//#endif //USE_SOLVER_FCBIILU2
