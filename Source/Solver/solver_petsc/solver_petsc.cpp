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
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    return (isInitialized == PETSC_TRUE);
}

void SolverInitializePetsc(int *argc, char ***argv, const char *file_options) 
{
    if (!SolverIsInitializedPetsc()) 
    {
        PetscErrorCode ierr = PetscInitialize(argc, argv, file_options, "solver");
        //prevent petsc from catching signals, petsc error handling is misleading for
        //unexperienced user and results in the opinion that errors emerge from petsc
        PetscPopSignalHandler();
        if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    }
}

bool SolverIsFinalizedPetsc() 
{
    PetscBool isFinalized;
    PetscErrorCode ierr = PetscFinalized(&isFinalized);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    return (isFinalized == PETSC_TRUE);
}

void SolverFinalizePetsc() 
{
    PetscErrorCode ierr = PetscFinalize();
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}

void MatrixDestroyDataPetsc(Mat **matrix) 
{
    PetscErrorCode ierr = MatDestroy(*matrix);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    delete *matrix;
    *matrix = NULL;
}


void MatrixInitDataPetsc(Mat **matrix, INMOST_MPI_Comm comm, const char *name) 
{
    PetscErrorCode ierr;
    if (matrix == NULL) throw INMOST::DataCorruptedInSolver;
    if (*matrix != NULL)
        MatrixDestroyDataPetsc(matrix);
    *matrix = new Mat();
#if !defined(USE_MPI)
    ierr = MatCreate(PETSC_COMM_WORLD, *matrix);
#else
    ierr = MatCreate(comm, *matrix);
#endif
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    ierr = MatSetOptionsPrefix(**matrix, name);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    ierr = MatSetFromOptions(**matrix);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}

void MatrixCopyDataPetsc(Mat **matrix, Mat *other_matrix) 
{
    PetscErrorCode ierr;
    if (matrix == NULL || other_matrix == NULL) throw INMOST::DataCorruptedInSolver;
    if (*matrix != NULL)
        MatrixDestroyDataPetsc(matrix);
    *matrix = new Mat();
    ierr = MatConvert(*other_matrix, MATSAME, MAT_INITIAL_MATRIX, *matrix);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}

void MatrixAssignDataPetsc(Mat *matrix, Mat *other_matrix) 
{
    PetscErrorCode ierr;
    ierr = MatCopy(*other_matrix, *matrix, DIFFERENT_NONZERO_PATTERN);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}


void
MatrixPreallocatePetsc(Mat *matrix, int local_size, int global_size, int *diag_nonzeroes, int *off_diag_nonzeroes) 
{
    PetscErrorCode ierr;
    ierr = MatSetSizes(*matrix, local_size, local_size, global_size, global_size);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    ierr = MatXAIJSetPreallocation(*matrix, 1, diag_nonzeroes, off_diag_nonzeroes, PETSC_NULL, PETSC_NULL);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}

void MatrixFillPetsc(Mat *matrix, int row, int cols, int *col_positions, double *col_values)
{
    PetscErrorCode ierr;
    ierr = MatSetValues(*matrix, 1, &row, cols, col_positions, col_values, INSERT_VALUES);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}

void MatrixFinalizePetsc(Mat *matrix) 
{
    PetscErrorCode ierr;
    ierr = MatAssemblyBegin(*matrix, MAT_FINAL_ASSEMBLY);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    ierr = MatAssemblyEnd(*matrix, MAT_FINAL_ASSEMBLY);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}


void VectorDestroyDataPetsc(Vec **vector) 
{
    PetscErrorCode ierr = VecDestroy(*vector);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    delete *vector;
    *vector = NULL;
}


void VectorInitDataPetsc(Vec **vector, INMOST_MPI_Comm comm, const char *name) 
{
    PetscErrorCode ierr;
    if (vector == NULL) throw INMOST::DataCorruptedInSolver;
    if (*vector != NULL)
        VectorDestroyDataPetsc(vector);
    *vector = new Vec();
#if !defined(USE_MPI)
    ierr = VecCreate(PETSC_COMM_WORLD, *vector);
#else
    ierr = VecCreate(comm, *vector);
#endif
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    ierr = VecSetOptionsPrefix(**vector, name);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    ierr = VecSetFromOptions(**vector);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}

void VectorCopyDataPetsc(Vec **vector, Vec *other_vector) 
{
    PetscErrorCode ierr;
    if (vector == NULL || other_vector == NULL) throw INMOST::DataCorruptedInSolver;
    if (*vector != NULL)
        VectorDestroyDataPetsc(vector);
    *vector = new Vec();
    ierr = VecDuplicate(*other_vector, *vector);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    ierr = VecCopy(*other_vector, **vector);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}


void VectorAssignDataPetsc(Vec *vector, Vec *other_vector) 
{
    PetscErrorCode ierr;
    if (vector == NULL || other_vector == NULL) throw INMOST::DataCorruptedInSolver;
    ierr = VecCopy(*other_vector, *vector);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}

void VectorPreallocatePetsc(Vec *vector, int local_size, int global_size) 
{
    PetscErrorCode ierr = VecSetSizes(*vector, local_size, global_size);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}

void VectorFillPetsc(Vec *vector, int size, int *positions, double *values) 
{
    PetscErrorCode ierr = VecSetValues(*vector, size, positions, values, INSERT_VALUES);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}

void VectorLoadPetsc(Vec *vector, int size, int *positions, double *values) 
{
    PetscErrorCode ierr = VecGetValues(*vector, size, positions, values);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}

void VectorFinalizePetsc(Vec *vector) 
{
    PetscErrorCode ierr;
    ierr = VecAssemblyBegin(*vector);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    ierr = VecAssemblyEnd(*vector);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}

void SolverDestroyDataPetsc(KSP **ksp) 
{
    PetscErrorCode ierr = KSPDestroy(*ksp);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    delete *ksp;
    *ksp = NULL;
}


void SolverInitDataPetsc(KSP **ksp, INMOST_MPI_Comm comm, const char *name) 
{
    PetscErrorCode ierr;
    if (ksp == NULL) throw INMOST::DataCorruptedInSolver;
    if (*ksp != NULL)
        SolverDestroyDataPetsc(ksp);
    *ksp = new KSP();
#if !defined(USE_MPI)
    ierr = KSPCreate(PETSC_COMM_WORLD, *ksp);
#else
    ierr = KSPCreate(comm, *ksp);
#endif
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    ierr = KSPSetOptionsPrefix(**ksp, name);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    ierr = KSPSetFromOptions(**ksp);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}

void SolverCopyDataPetsc(KSP **ksp, KSP *other_ksp, INMOST_MPI_Comm comm) 
{
    PetscErrorCode ierr;
    if (ksp == NULL) throw INMOST::DataCorruptedInSolver;
    if (*ksp != NULL)
        SolverDestroyDataPetsc(ksp);
    *ksp = new KSP();
#if !defined(USE_MPI)
    ierr = KSPCreate(PETSC_COMM_WORLD, *ksp);
#else
    ierr = KSPCreate(comm, *ksp);
#endif
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    char *prefix;
    ierr = KSPGetOptionsPrefix(*other_ksp, const_cast<const char **>(&prefix));
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    ierr = KSPSetOptionsPrefix(**ksp, prefix);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    ierr = KSPSetFromOptions(**ksp);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}

void SolverAssignDataPetsc(KSP *ksp, KSP *other_ksp)
{
    PetscErrorCode ierr;
    char *prefix;
    ierr = KSPGetOptionsPrefix(*other_ksp, const_cast<const char **>(&prefix));
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    ierr = KSPSetOptionsPrefix(*ksp, prefix);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    ierr = KSPSetFromOptions(*ksp);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}


void SolverSetMatrixPetsc(KSP *ksp, Mat *matrix, bool same_pattern, bool reuse_preconditioner) 
{
    PetscErrorCode ierr;
    if (reuse_preconditioner)
        ierr = KSPSetReusePreconditioner(*ksp, PETSC_TRUE);
    else
        ierr = KSPSetReusePreconditioner(*ksp, PETSC_FALSE);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    ierr = KSPSetOperators(*ksp, *matrix,
                           *matrix);//,reuse_preconditioner? SAME_PRECONDITIONER : (same_pattern? SAME_NONZERO_PATTERN : DIFFERENT_NONZERO_PATTERN));
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}

bool SolverSolvePetsc(KSP *ksp, Vec *rhs, Vec *sol) 
{
    PetscErrorCode ierr;
    PetscBool haveksp;
    bool guess = true;
    char typeksp[2048];
    char *prefix;
    ierr = KSPGetOptionsPrefix(*ksp, const_cast<const char **>(&prefix));
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
//due to https://www.mcs.anl.gov/petsc/documentation/changes/37.html
#if PETSC_VERSION_MAJOR > 3 || (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 7)
    ierr = PetscOptionsGetString(NULL,prefix, "-ksp_type", typeksp, 2048, &haveksp);
#else //PETSC_VERSION
    ierr = PetscOptionsGetString(prefix, "-ksp_type", typeksp, 2048, &haveksp);
#endif //PETSC_VERSION
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    if (haveksp && !strcmp(typeksp, "preonly")) guess = false;
    if (guess) 
    {
        ierr = KSPSetInitialGuessNonzero(*ksp, PETSC_TRUE);
        if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    }
    //MatNullSpace nullsp;
    //ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL,&nullsp);
    //ierr = KSPSetNullSpace(*ksp,nullsp);
    ierr = KSPSolve(*ksp, *rhs, *sol);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(*ksp, &reason);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    return reason >= 0;
}

int SolverIterationNumberPetsc(KSP *ksp) 
{
    PetscInt its;
    PetscErrorCode ierr = KSPGetIterationNumber(*ksp, &its);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    return its;
}

double SolverResidualNormPetsc(KSP *ksp) 
{
    PetscReal norm;
    PetscErrorCode ierr = KSPGetResidualNorm(*ksp, &norm);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    return norm;
}


void SolverSetTolerancesPetsc(KSP *ksp, double rtol, double atol, double divtol, int maxits) 
{
    PetscErrorCode ierr = KSPSetTolerances(*ksp, rtol, atol, divtol, maxits);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}


void SolverSetOverlapPetsc(KSP *ksp, int levels) 
{
    PetscErrorCode ierr;
    PC pc;
    PCType pc_type;
    ierr = KSPGetPC(*ksp, &pc);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    ierr = PCGetType(pc, &pc_type);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    if (!strcmp(pc_type, "asm")) 
    {
        ierr = PCASMSetOverlap(pc, levels);
    } 
    else if (!strcmp(pc_type, "gasm")) 
    {
        ierr = PCGASMSetOverlap(pc, levels);
    }
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}

void SolverSetDropTolerancePetsc(KSP *ksp, double dtol) 
{
    PetscErrorCode ierr;
    PC pc;
    ierr = KSPGetPC(*ksp, &pc);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    ierr = PCFactorSetDropTolerance(pc, dtol, 0.01, 10000); //2 other parameters are set to extreme
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}

void SolverSetFillLevelPetsc(KSP *ksp, double lfill) 
{
    PetscErrorCode ierr;
    PC pc;
    ierr = KSPGetPC(*ksp, &pc);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    ierr = PCFactorSetLevels(pc, static_cast<PetscInt>(lfill));
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
}

const char *SolverConvergedReasonPetsc(KSP *ksp) 
{
    static char reason_str[4096];
    KSPConvergedReason reason;
    PetscErrorCode ierr = KSPGetConvergedReason(*ksp, &reason);
    if (ierr != PETSC_SUCCESS) throw INMOST::ErrorInSolver;
    switch (reason) 
    {
        case KSP_CONVERGED_RTOL:
        case KSP_CONVERGED_RTOL_NORMAL:
            strcpy(reason_str, "norm decreased by a factor of relative tolerance");
            break;
        case KSP_CONVERGED_ATOL:
        case KSP_CONVERGED_ATOL_NORMAL:
            strcpy(reason_str, "norm less then absolute tolerance");
            break;
        case KSP_CONVERGED_ITS:
            strcpy(reason_str, "converged by direct solver");
            break;
        case KSP_CONVERGED_CG_NEG_CURVE:
        case KSP_CONVERGED_CG_CONSTRAINED:
            strcpy(reason_str, "converged due to some condition in conjugate gradient method");
            break;
        case KSP_CONVERGED_STEP_LENGTH:
            strcpy(reason_str, "converged due to step length");
            break;
        case KSP_CONVERGED_HAPPY_BREAKDOWN :
            strcpy(reason_str, "converged due to happy breakdown");
            break;
        case KSP_CONVERGED_ITERATING:
            strcpy(reason_str, "converged");
            break;
        case KSP_DIVERGED_NULL:
            strcpy(reason_str, "diverged due to null");
            break;
        case KSP_DIVERGED_ITS:
            strcpy(reason_str, "diverged due to maximum iteration number");
            break;
        case KSP_DIVERGED_DTOL:
            strcpy(reason_str, "diverged as divergence tolerance was reached");
            break;
        case KSP_DIVERGED_BREAKDOWN:
            strcpy(reason_str, "diverged due to breakdown in method");
            break;
        case KSP_DIVERGED_BREAKDOWN_BICG:
            strcpy(reason_str, "diverged due to breakdown in Biconjugate Gradients method");
            break;
        case KSP_DIVERGED_NONSYMMETRIC:
            strcpy(reason_str, "diverged since matrix is nonsymmetric");
            break;
        case KSP_DIVERGED_INDEFINITE_PC:
            strcpy(reason_str, "diverged since preconditioner is not defined");
            break;
        case KSP_DIVERGED_NANORINF:
            strcpy(reason_str, "diverged since not a number or infinite was encountered");
            break;
        case KSP_DIVERGED_INDEFINITE_MAT:
            strcpy(reason_str, "diverged due to matrix is not definite");
            break;
        default:
            strcpy(reason_str, "reason code was not defined in manual");
            break;
    }
    return reason_str;
}

#endif
#endif


