#include "solver_ani.h"

#if defined(_WIN32) || defined(_WIN64) || defined(WIN32) || defined(WIN64)
__declspec( dllexport ) void solver_ani_stub(){} //to avoid LNK4221 warning
#endif

#if defined(USE_SOLVER_ANI)

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

#if defined (__APPLE__) || defined(MAXOSX) || defined(__linux__)
#define ILUOO    iluoo_
#define DDOT    ddot_
#define ILU0    ilu0_
#define PREVEC0    prevec0_
#define PREVEC2    prevec2_
#define MATVEC    matvec_
#define SLPBCGS    slpbcgs_
#endif
#if defined(_WIN32) || defined(_WIN64) || defined(WIN32) || defined(WIN64)
#define ILUOO	ILUOO
#define DDOT	DDOT
#define ILU0	ILU0
#define PREVEC0	PREVEC0
#define PREVEC2	PREVEC2
#define MATVEC	MATVEC
#define SLPBCGS	SLPBCGS
#endif


extern "C"
{
//	 from blas 
double DDOT(int *n, double *dx, int *incx, double *dy, int *incy);
//	 from aniILU 
void ILU0(int *n, double *a, int *ja, int *ia, double *alu, int *jlu, int *ju, int *iw, int *ierr);
void ILUOO(int *n, int *ia, int *ja, double *a, double *tau1, double *tau2, int *verb,
           double *rW, int *iW, int *MaxWr, int *MaxWi, double *partlur, double *partlurout,
           int *UsedWr, int *UsedWi, int *ierr);
void PREVEC0(int *iprevec, int *ichange, double *x, double *y, int *iW, double *rW);
void PREVEC2(int *iprevec, int *ichange, double *x, double *y, int *iW, double *rW);
void MATVEC(int *imatvec, double *alpha, double *x, double *beta, double *y, int *ia, int *ja, double *a);
void SLPBCGS(void (*prevec)(int *, int *, double *, double *, int *, double *), int *iprevec, int *iW, double *rW,
             void (*matvec)(int *, double *, double *, double *, double *, int *, int *, double *),
             int *imatvec, int *ia, int *ja, double *a,
             double *work, int *mw, int *nw,
             int *n, double *rhs, double *sol,
             int *iter, double *resid,
             int *info, int *nuint);
}

int ilu2 = 1;
double set_epsilon = 1e-5;

/* for internal use */
static double residual_fortran_indices(matrix *A, vector *x, vector *b) {
    double d = 0.0, c;
    int i, j;
    for (i = 0; i < A->n; i++, d += c * c)
        for (c = b->v[i], j = A->ia[i] - 1; j < A->ia[i + 1] - 1; j++)
            c -= A->A[j] * x->v[A->ja[j] - 1];
    return sqrt(d);
}


/* Initialize solver with new matrix A */
int newmatrixbcg(bcg *s, matrix *A, bool same_precond) {
    if (s->n != A->n && same_precond) throw INMOST::CannotReusePreconditionerOfDifferentSize;
    s->n = A->n;
    s->nz = A->ia[A->n] - A->ia[0];
    s->ia = A->ia, s->ja = A->ja;
    if (!same_precond) {
        if (!ilu2) {
            s->rW = (double *) realloc(s->rW, sizeof(double) * 2 * (s->nz + s->n * 8));
//			memset(s->rW,0,sizeof(double)*2*(s->nz + s->n*8    ));
            s->iW = (int *) realloc(s->iW, sizeof(int) * 2 * (s->nz + s->n * 2 + 1));
//			memset(s->iW,0,sizeof(int)   *2*(s->nz + s->n*2 + 1));
            s->ipBCG = s->nz;
        } else {
            s->rW = (double *) realloc(s->rW, sizeof(double) * 2 * (s->nz * 20)); // ilu2
//			memset(s->rW,0,sizeof(double)*2*(s->nz*20));
            s->iW = (int *) realloc(s->iW, sizeof(int) * 2 * (s->nz * 25)); // ilu2
//			memset(s->iW,0,sizeof(int)   *2*(s->nz*25));
        }
        s->ipjLU = s->n + 1;
        s->ipiw = s->ipjLU + s->nz;
        return renewbcg(s, A->A);
    } else return 0;
}

/* solver */
/* Initialize bcg solver */
int initbcg(bcg *s, matrix *A) {
    s->eps = set_epsilon;
    s->rW = NULL;
    s->iW = NULL;
    return newmatrixbcg(s, A, false);
}

/* Reinitialize solver preconditioner with new matrix A */
int renewbcg(bcg *s, double *A) {
    int ierr = 0;
    int v, k;

    /* Oh, we need to check the indexing again */
    v = 1 - s->ia[0];
    if (v) for (k = 0; k <= s->n; k++) s->ia[k] += v;
    if (v) for (k = 0; k < s->nz; k++) s->ja[k] += v;
    /* Reinitialize ilu0 */
    s->A = A;

    if (!ilu2)
        ILU0(&s->n, s->A, s->ja, s->ia, s->rW, s->iW + s->ipjLU, s->iW, s->iW + s->ipiw, &ierr);
    else {
        int verb = 0, MaxWr = 2 * s->nz * 20, MaxWi = 2 * s->nz * 25, UsedWr, UsedWi;
        double tau1 = 1e-3, tau2 = 1e-6, partlur = 0.5, partlurout;
        ILUOO(&s->n, s->ia, s->ja, s->A, &tau1, &tau2, &verb, s->rW, s->iW,
              &MaxWr, &MaxWi, &partlur, &partlurout, &UsedWr, &UsedWi, &ierr);
        s->ipBCG = UsedWr;
    }
    if (ierr) printf("initialization of ilu failed, zero pivot=%d\n", ierr);
    return ierr;
}

/* Solve linear system */
int solvebcg_ani(bcg *s, vector *b, vector *x) {
    int NUNIT = 6, INFO = 0, eight = 8;

    s->ITER = 1000;
    matrix temp;
    temp.n = s->n;
    temp.ia = s->ia;
    temp.ja = s->ja;
    temp.A = s->A;
    s->RESID = s->eps * residual_fortran_indices(&temp, x, b);
    if (!ilu2)
        SLPBCGS(PREVEC0, &s->n, s->iW, s->rW,
                MATVEC, &s->n, s->ia, s->ja, s->A,
                s->rW + s->ipBCG, &s->n, &eight,
                &s->n, b->v, x->v,
                &s->ITER, &s->RESID,
                &INFO, &NUNIT);
    else
        SLPBCGS(PREVEC2, &s->n, s->iW, s->rW,
                MATVEC, &s->n, s->ia, s->ja, s->A,
                s->rW + s->ipBCG, &s->n, &eight,
                &s->n, b->v, x->v,
                &s->ITER, &s->RESID,
                &INFO, &NUNIT);
    return INFO;
}

/* Free memory used by solver */
void freebcg(bcg *s) {
    free(s->rW), free(s->iW);
}

#endif //USE_SOLVER_ANI
