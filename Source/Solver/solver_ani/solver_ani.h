#ifndef SOLVER_ANI_H_INCLUDED
#define SOLVER_ANI_H_INCLUDED
#include "inmost_solver.h"

/* BiCGStab solver structure */
typedef struct
{
    int n, nz, ierr, ipjLU, ipiw, ipBCG;
    double RESID;
    int ITER, INFO, NUNIT;
    int *iW;
    double *rW;
    int *ia, *ja;
    double *A;
    double eps;
} bcg;

typedef struct
{
    int n;
    int * ia;
    int * ja;
    double * A;
} matrix;

typedef struct
{
    int n;
    double * v;
} vector;

/* Initialize bcg solver */
int initbcg(bcg *s, matrix *A);

/* Reinitialize solver preconditioner with new matrix A */
int renewbcg(bcg *s, double *A);

/* Solve linear system */
int solvebcg_ani(bcg *s, vector *b, vector *x);

/* Free memory used by solver */
void freebcg(bcg *s);

/* Initialize solver with new matrix A */
int newmatrixbcg(bcg *s, matrix *A, bool same_precond);


#endif //SOLVER_ANI_H_INCLUDED
