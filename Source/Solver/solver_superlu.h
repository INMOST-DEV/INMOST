#ifndef _SOLVER_SUPERLU
#define _SOLVER_SUPERLU

//IF CERTAIN FUNCTIONALITY IS NOT AVAILABLE, PLEASE THROW INMOST::NotImplemented EXCEPTION

void SolverInitializeSuperLU(int * argc,char *** argv, const char * file_options);
void MatrixFillSuperLU(void * data, int size, int nnz, int * col_ranges, int * col_positions, double * col_values);
void SolverDestroyDataSuperLU(void ** data);
void SolverInitDataSuperLU(void ** data);
bool SolverSolveSuperLU(void * data, void * rhs_data, void * sol_data);
const char * SolverConvergedReasonSuperLU(void * data);

#endif
