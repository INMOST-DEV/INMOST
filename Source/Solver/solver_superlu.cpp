#define _CRT_SECURE_NO_WARNINGS
#include "inmost_solver.h"
#if defined(USE_SOLVER)
#if defined(USE_SOLVER_SUPERLU)
#include "superlu/slu_ddefs.h"
#include "solver_superlu.h"

struct SuperLU_data
{
	SuperMatrix A, L, U;
	int * perm_r;
	int * perm_c;
	int size, info;
	superlu_options_t options;
	SuperLUStat_t stat;
};

void SolverInitializeSuperLU(int * argc,char *** argv, const char * file_options)
{
	//read options from file and arguments
}


void MatrixFillSuperLU(void * data, int size, int nnz, int * col_ranges, int * col_positions, double * col_values)
{
	SuperLU_data * m = static_cast<SuperLU_data *>(data);
	dCreate_CompCol_Matrix(&m->A,size,size,nnz,col_values,col_positions,col_ranges,SLU_NR,SLU_D,SLU_GE);
	m->size = size;
	m->perm_c = new int [size];
	m->perm_r = new int [size];
}


void SolverDestroyDataSuperLU(void ** data)
{
	SuperLU_data * m = static_cast<SuperLU_data *>(*data);
	Destroy_CompCol_Matrix(&m->A);
	Destroy_SuperNode_Matrix(&m->L);
	Destroy_CompCol_Matrix(&m->U);
	StatFree(&m->stat);
	if( m->perm_c != NULL ) delete [] m->perm_c;
	if( m->perm_r != NULL ) delete [] m->perm_r;
	delete m;
	*data = NULL;
}

void SolverInitDataSuperLU(void ** data)
{
	if( data == NULL ) throw INMOST::DataCorruptedInSolver;
	if( *data != NULL ) 
		SolverDestroyDataSuperLU(data);
	*data = static_cast<void *>(new SuperLU_data);
	SuperLU_data * m = static_cast<SuperLU_data *>(*data);
	set_default_options(&(m->options));
	StatInit(&(m->stat));
	m->perm_c = NULL;
	m->perm_r = NULL;
	m->size = 0;
}


bool SolverSolveSuperLU(void * data, void * rhs_data, void * sol_data)
{
	SuperLU_data * m = static_cast<SuperLU_data *>(data);
	assert(m->size != 0);
	SuperMatrix B;
	memcpy(sol_data,rhs_data,sizeof(double)*m->size);
	dCreate_Dense_Matrix(&B,m->size,1,(double *)sol_data,m->size,SLU_DN,SLU_D,SLU_GE);
	dgssv(&m->options,&m->A,m->perm_c,m->perm_r,&m->L,&m->U,&B,&m->stat,&m->info);
	Destroy_SuperMatrix_Store(&B);
	return m->info == 0;
}

const char * SolverConvergedReasonSuperLU(void * data)
{
	SuperLU_data * m = static_cast<SuperLU_data *>(data);
	static char reason_str[4096];
	if( m->info <= m->size )
		sprintf(reason_str,"diagonal entry of U-factor is exactly singular at %d/%d",m->info,m->size);
	else if( m->info > m->size )
		sprintf(reason_str,"memory allocation failed after %d bytes were allocated",m->info-m->size);
	else if( m->info == 0 )
		strcpy(reason_str,"factorization was successfull");
	else
		sprintf(reason_str,"unknown exit code %d",m->info);
	return reason_str;
}

#endif
#endif
