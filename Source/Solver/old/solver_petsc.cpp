#define _CRT_SECURE_NO_WARNINGS
#include "inmost_solver.h"
#if defined(USE_SOLVER)
#if defined(USE_SOLVER_PETSC)
#include <petsc.h>
#include "solver_petsc.h"
#define PETSC_SUCCESS 0


bool SolverIsInitializedPetsc()
{
	PetscBool isInitialized;
	PetscErrorCode ierr = PetscInitialized(&isInitialized);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	return (isInitialized == PETSC_TRUE);
}

void SolverInitializePetsc(int * argc,char *** argv, const char * file_options)
{
	if( !SolverIsInitializedPetsc() )
	{
		PetscErrorCode ierr = PetscInitialize(argc,argv,file_options,"solver");
		//prevent petsc from catching signals, petsc error handling is misleading for 
		//unexperienced user and results in the opinion that errors emerge from petsc
		PetscPopSignalHandler(); 
		if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	}
}

bool SolverIsFinalizedPetsc()
{
	PetscBool isFinalized;
	PetscErrorCode ierr = PetscFinalized(&isFinalized);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	return (isFinalized == PETSC_TRUE);
}

void SolverFinalizePetsc() {
    if (!SolverIsFinalizedPetsc()) {
        PetscErrorCode ierr = PetscFinalize();
        if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    }
}

void MatrixDestroyDataPetsc(void ** data)
{
	Mat * m = static_cast<Mat *>(*data);
	PetscErrorCode ierr = MatDestroy(m);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	delete m;
	*data = NULL;
}


void MatrixInitDataPetsc(void ** data, INMOST_MPI_Comm comm, const char * name)
{
	PetscErrorCode ierr;
	if( data == NULL ) throw INMOST::DataCorruptedInSolver;
	if( *data != NULL )
		MatrixDestroyDataPetsc(data);
	*data = static_cast<void *>(new Mat);
#if !defined(USE_MPI)
	ierr = MatCreate(PETSC_COMM_WORLD,static_cast<Mat *>(*data));
#else
	ierr = MatCreate(comm,static_cast<Mat *>(*data));
#endif
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	ierr = MatSetOptionsPrefix(*static_cast<Mat *>(*data),name);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	ierr = MatSetFromOptions(*static_cast<Mat *>(*data));
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}

void MatrixCopyDataPetsc(void ** data, void * other_data)
{
	PetscErrorCode ierr;
	if( data == NULL || other_data == NULL ) throw INMOST::DataCorruptedInSolver;
	if( *data != NULL )
		MatrixDestroyDataPetsc(data);
	*data = static_cast<void *>(new Mat);
	ierr = MatConvert(*static_cast<Mat *>(other_data),MATSAME,MAT_INITIAL_MATRIX,static_cast<Mat *>(*data));
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}

void MatrixAssignDataPetsc(void * data, void * other_data)
{
	PetscErrorCode ierr;
	if( data == NULL || other_data == NULL ) throw INMOST::DataCorruptedInSolver;
	ierr = MatCopy(*static_cast<Mat *>(other_data),*static_cast<Mat *>(data),DIFFERENT_NONZERO_PATTERN);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}


void MatrixPreallocatePetsc(void * data, int local_size, int global_size, int * diag_nonzeroes, int * off_diag_nonzeroes)
{
	PetscErrorCode ierr;
	Mat * m = static_cast<Mat *>(data);
	ierr = MatSetSizes(*m,local_size,local_size,global_size,global_size);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	ierr = MatXAIJSetPreallocation(*m,1,diag_nonzeroes,off_diag_nonzeroes,PETSC_NULL,PETSC_NULL);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}

void MatrixFillPetsc(void * data, int row, int cols, int * col_positions, double * col_values)
{
	PetscErrorCode ierr;
	Mat * m = static_cast<Mat *>(data);
	ierr = MatSetValues(*m,1,&row,cols,col_positions,col_values,INSERT_VALUES);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}

void MatrixFinalizePetsc(void * data)
{
	PetscErrorCode ierr;
	Mat * m = static_cast<Mat *>(data);
	ierr = MatAssemblyBegin(*m,MAT_FINAL_ASSEMBLY);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	ierr = MatAssemblyEnd(*m,MAT_FINAL_ASSEMBLY);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}


void VectorDestroyDataPetsc(void ** data)
{
	Vec * m = static_cast<Vec *>(*data);
	PetscErrorCode ierr = VecDestroy(m);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	delete m;
	*data = NULL;
}


void VectorInitDataPetsc(void ** data, INMOST_MPI_Comm comm, const char * name)
{
	PetscErrorCode ierr;
	if( data == NULL ) throw INMOST::DataCorruptedInSolver;
	if( *data != NULL ) 
		VectorDestroyDataPetsc(data);
	*data = static_cast<void *>(new Vec);
#if !defined(USE_MPI)
	ierr = VecCreate(PETSC_COMM_WORLD,static_cast<Vec *>(*data));
#else
	ierr = VecCreate(comm,static_cast<Vec *>(*data));
#endif
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	ierr = VecSetOptionsPrefix(*static_cast<Vec *>(*data),name);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	ierr = VecSetFromOptions(*static_cast<Vec *>(*data));
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}

void VectorCopyDataPetsc(void ** data, void * other_data)
{
	PetscErrorCode ierr;
	if( data == NULL || other_data == NULL ) throw INMOST::DataCorruptedInSolver;
	if( *data != NULL )
		VectorDestroyDataPetsc(data);
	*data = static_cast<void *>(new Vec);
	ierr = VecDuplicate(*static_cast<Vec *>(other_data),static_cast<Vec *>(*data));
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	ierr = VecCopy(*static_cast<Vec *>(other_data),*static_cast<Vec *>(*data));
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}

void VectorAssignDataPetsc(void * data, void * other_data)
{
	PetscErrorCode ierr;
	if( data == NULL || other_data == NULL ) throw INMOST::DataCorruptedInSolver;
	ierr = VecCopy(*static_cast<Vec *>(other_data),*static_cast<Vec *>(data));
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}



void VectorPreallocatePetsc(void * data, int local_size, int global_size)
{
	Vec * m = static_cast<Vec *>(data);
	PetscErrorCode ierr = VecSetSizes(*m,local_size,global_size);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}
void VectorFillPetsc(void * data, int size, int * positions, double * values)
{
	Vec * m = static_cast<Vec *>(data);
	PetscErrorCode ierr = VecSetValues(*m,size,positions,values,INSERT_VALUES);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}
void VectorLoadPetsc(void * data, int size, int * positions, double * values)
{
	Vec * m = static_cast<Vec *>(data);
	PetscErrorCode ierr = VecGetValues(*m,size,positions,values);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}

void VectorFinalizePetsc(void * data)
{
	PetscErrorCode ierr;
	Vec * m = static_cast<Vec *>(data);
	ierr = VecAssemblyBegin(*m);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	ierr = VecAssemblyEnd(*m);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}

void SolverDestroyDataPetsc(void ** data)
{
	KSP * m = static_cast<KSP *>(*data);
	PetscErrorCode ierr = KSPDestroy(m);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	delete m;
	*data = NULL;
}


void SolverInitDataPetsc(void ** data, INMOST_MPI_Comm comm, const char * name)
{
	PetscErrorCode ierr;
	if( data == NULL ) throw INMOST::DataCorruptedInSolver;
	if( *data != NULL )
		SolverDestroyDataPetsc(data);
	*data = static_cast<void *>(new KSP);
#if !defined(USE_MPI)
	ierr = KSPCreate(PETSC_COMM_WORLD,static_cast<KSP *>(*data));
#else
	ierr = KSPCreate(comm,static_cast<KSP *>(*data));
#endif
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	ierr = KSPSetOptionsPrefix(*static_cast<KSP *>(*data),name);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	ierr = KSPSetFromOptions(*static_cast<KSP *>(*data));
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}

void SolverCopyDataPetsc(void ** data, void * other_data, INMOST_MPI_Comm comm)
{
	PetscErrorCode ierr;
	if( data == NULL || other_data == NULL ) throw INMOST::DataCorruptedInSolver;
	if( *data != NULL )
		SolverDestroyDataPetsc(data);
	*data = static_cast<void *>(new KSP);
#if !defined(USE_MPI)
	ierr = KSPCreate(PETSC_COMM_WORLD,static_cast<KSP *>(*data));
#else
	ierr = KSPCreate(comm,static_cast<KSP *>(*data));
#endif
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	KSP * mine = static_cast<KSP *>(*data);
	KSP * other = static_cast<KSP *>(other_data);
	char * prefix;
	ierr = KSPGetOptionsPrefix(*other,const_cast<const char **>(&prefix));
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	ierr = KSPSetOptionsPrefix(*mine,prefix);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	ierr = KSPSetFromOptions(*mine);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}
void SolverAssignDataPetsc(void * data, void * other_data)
{
	PetscErrorCode ierr;
	if( data == NULL || other_data == NULL ) throw INMOST::DataCorruptedInSolver;
	KSP * mine = static_cast<KSP *>(data);
	KSP * other = static_cast<KSP *>(other_data);
	char * prefix;
	ierr = KSPGetOptionsPrefix(*other,const_cast<const char **>(&prefix));
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	ierr = KSPSetOptionsPrefix(*mine,prefix);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	ierr = KSPSetFromOptions(*mine);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}


void SolverSetMatrixPetsc(void * data, void * matrix_data, bool same_pattern, bool reuse_preconditioner)
{
	KSP * m = static_cast<KSP *>(data);
	Mat * mat = static_cast<Mat *>(matrix_data);
	PetscErrorCode ierr;
	if( reuse_preconditioner ) 
		ierr = KSPSetReusePreconditioner(*m,PETSC_TRUE);
	else 
		ierr = KSPSetReusePreconditioner(*m,PETSC_FALSE);
	ierr = KSPSetOperators(*m,*mat,*mat);//,reuse_preconditioner? SAME_PRECONDITIONER : (same_pattern? SAME_NONZERO_PATTERN : DIFFERENT_NONZERO_PATTERN));
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}

bool SolverSolvePetsc(void * data, void * rhs_data, void * sol_data)
{
	PetscErrorCode ierr;
	PetscBool haveksp;
	Vec * rhs = static_cast<Vec *>(rhs_data), * sol = static_cast<Vec *>(sol_data);
	KSP * ksp = static_cast<KSP *>(data);
	bool guess = true;
	char typeksp[2048];
	char * prefix;
	ierr = KSPGetOptionsPrefix(*ksp,const_cast<const char **>(&prefix));
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	ierr = PetscOptionsGetString(prefix,"-ksp_type",typeksp,2048,&haveksp);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	if( haveksp && !strcmp(typeksp,"preonly") ) guess = false;
	if( guess )
	{
		ierr = KSPSetInitialGuessNonzero(*ksp,PETSC_TRUE);
		if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	}
  //MatNullSpace nullsp;
  //ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL,&nullsp);
  //ierr = KSPSetNullSpace(*ksp,nullsp);
	ierr = KSPSolve(*ksp,*rhs,*sol);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	KSPConvergedReason reason;
	ierr = KSPGetConvergedReason(*ksp,&reason);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	if( reason >= 0 )
		return true;
	return false;
}
int SolverIterationNumberPetsc(void * data)
{
	PetscInt its;
	KSP * ksp = static_cast<KSP *>(data);
	PetscErrorCode ierr = KSPGetIterationNumber(*ksp,&its);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	return its;
}
double SolverResidualNormPetsc(void * data)
{
	PetscReal norm;
	KSP * ksp = static_cast<KSP *>(data);
	PetscErrorCode ierr = KSPGetResidualNorm(*ksp,&norm);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	return norm;
}


void SolverSetTolerancesPetsc(void * data, double rtol, double atol, double divtol, int maxits)
{
	KSP * ksp = static_cast<KSP *>(data);
	PetscErrorCode ierr = KSPSetTolerances(*ksp,rtol,atol,divtol,maxits);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}


void SolverSetOverlapPetsc(void * data, int levels)
{
	KSP * ksp = static_cast<KSP *>(data);
	PetscErrorCode ierr;
	PC pc; PCType pc_type;
	ierr = KSPGetPC(*ksp,&pc);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	ierr = PCGetType(pc,&pc_type);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	if( !strcmp(pc_type,"asm") )
	{
		ierr = PCASMSetOverlap(pc,levels);
	}
	else if( !strcmp(pc_type,"gasm") )
	{
		ierr = PCGASMSetOverlap(pc,levels);
	}
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}

void SolverSetDropTolerancePetsc(void * data, double dtol)
{
	KSP * ksp = static_cast<KSP *>(data);
	PetscErrorCode ierr;
	PC pc;
	ierr = KSPGetPC(*ksp,&pc);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	ierr = PCFactorSetDropTolerance(pc,dtol,0.01,10000); //2 other parameters are set to extreme
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}

void SolverSetFillLevelPetsc(void * data, double lfill)
{
	KSP * ksp = static_cast<KSP *>(data);
	PetscErrorCode ierr;
	PC pc;
	ierr = KSPGetPC(*ksp,&pc);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	ierr = PCFactorSetLevels(pc,static_cast<PetscInt>(lfill));
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
}

const char * SolverConvergedReasonPetsc(void * data)
{
	static char reason_str[4096];
	KSPConvergedReason reason;
	KSP * ksp = static_cast<KSP *>(data);
	PetscErrorCode ierr = KSPGetConvergedReason(*ksp,&reason);
	if( ierr != PETSC_SUCCESS ) throw INMOST::ErrorInSolver;
	switch(reason)
	{
	case KSP_CONVERGED_RTOL:
	case KSP_CONVERGED_RTOL_NORMAL:
		strcpy(reason_str,"norm decreased by a factor of relative tolerance");
		break;
	case KSP_CONVERGED_ATOL:
	case KSP_CONVERGED_ATOL_NORMAL:
		strcpy(reason_str,"norm less then absolute tolerance");
		break;
	case  KSP_CONVERGED_ITS:
		strcpy(reason_str,"converged by direct solver");
		break;
	case KSP_CONVERGED_CG_NEG_CURVE:
	case KSP_CONVERGED_CG_CONSTRAINED:
		strcpy(reason_str,"converged due to some condition in conjugate gradient method");
		break;
	case KSP_CONVERGED_STEP_LENGTH:
		strcpy(reason_str,"converged due to step length");
		break;
	case KSP_CONVERGED_HAPPY_BREAKDOWN :
		strcpy(reason_str,"converged due to happy breakdown");
		break;
	case  KSP_CONVERGED_ITERATING:
		strcpy(reason_str,"converged");
		break;
	case KSP_DIVERGED_NULL:
		strcpy(reason_str,"diverged due to null");
		break;
	case KSP_DIVERGED_ITS:
		strcpy(reason_str,"diverged due to maximum iteration number");
		break;
	case KSP_DIVERGED_DTOL:
		strcpy(reason_str,"diverged as divergence tolerance was reached");
		break;
	case KSP_DIVERGED_BREAKDOWN:
		strcpy(reason_str,"diverged due to breakdown in method");
		break;
	case KSP_DIVERGED_BREAKDOWN_BICG:
		strcpy(reason_str,"diverged due to breakdown in Biconjugate Gradients method");
		break;
	case KSP_DIVERGED_NONSYMMETRIC:
		strcpy(reason_str,"diverged since matrix is nonsymmetric");
		break;
	case KSP_DIVERGED_INDEFINITE_PC:
		strcpy(reason_str,"diverged since preconditioner is not defined");
		break;
	case KSP_DIVERGED_NANORINF:
		strcpy(reason_str,"diverged since not a number or infinite was encountered");
		break;
	case KSP_DIVERGED_INDEFINITE_MAT:
		strcpy(reason_str,"diverged due to matrix is not definite");
		break;
	default:
		strcpy(reason_str,"reason code was not defined in manual");
		break;
	}
	return reason_str;
}
#endif
#endif
