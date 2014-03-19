#ifndef _SOLVER_PETSC
#define _SOLVER_PETSC

//IF CERTAIN FUNCTIONALITY IS NOT AVAILABLE, PLEASE THROW INMOST::NotImplemented EXCEPTION

void MatrixInitDataPetsc(void ** data, INMOST_MPI_Comm comm, const char * name);
void MatrixCopyDataPetsc(void ** data, void * other_data);
void MatrixAssignDataPetsc(void * data, void * other_data);
void MatrixDestroyDataPetsc(void ** data);
void MatrixPreallocatePetsc(void * data, int local_size, int global_size, int * diag_nonzeroes, int * off_diag_nonzeroes);
void MatrixFillPetsc(void * data, int row, int cols, int * col_positions, double * col_values);
void MatrixFinalizePetsc(void * data);


void VectorInitDataPetsc(void ** data, INMOST_MPI_Comm comm, const char * name);
void VectorCopyDataPetsc(void ** data, void * other_data);
void VectorAssignDataPetsc(void * data, void * other_data);
void VectorDestroyDataPetsc(void ** data);
void VectorPreallocatePetsc(void * data, int local_size, int global_size);
void VectorFillPetsc(void * data, int size, int * positions, double * values);
void VectorLoadPetsc(void * data, int size, int * positions, double * values);
void VectorFinalizePetsc(void * data);

bool SolverIsInitializedPetsc();
void SolverInitializePetsc(int * argc,char *** argv, const char * file_options);
bool SolverIsFinalizedPetsc();
void SolverFinalizePetsc();

void SolverInitDataPetsc(void ** data, INMOST_MPI_Comm comm, const char * name);
void SolverCopyDataPetsc(void ** data, void * other_data, INMOST_MPI_Comm comm);
void SolverAssignDataPetsc(void * data, void * other_data);
void SolverDestroyDataPetsc(void ** data);
void SolverSetMatrixPetsc(void * data, void * matrix_data, bool same_pattern, bool reuse_preconditioner);
bool SolverSolvePetsc(void * data, void * rhs_data, void * sol_data);
int SolverIterationNumberPetsc(void * data);
double SolverResidualNormPetsc(void * data);

#endif
