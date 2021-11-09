
#ifndef INMOST_COMMON_INCLUDED
#define INMOST_COMMON_INCLUDED


#include "inmost_options.h"
#if !defined(INMOST_OPTIONS_CMAKE_INCLUDED)
#pragma message("CMake options are not used")
//no options from cmake -- define minimum
//#define USE_OMP
//#define USE_OPENCL

#define USE_MESH

#define USE_PARTITIONER
//#define USE_PARTITIONER_ZOLTAN
//#define USE_PARTITIONER_PARMETIS
//#define USE_PARTITIONER_METIS

#define USE_SOLVER
//#define USE_SOLVER_METIS
//#define USE_SOLVER_PETSC
//#define USE_SOLVER_TRILINOS
//#define USE_SOLVER_ANI


#define USE_NONLINEAR
//#define USE_NONLINEAR_PETSC
//#define USE_NONLINEAR_TRILINOS
//#define USE_NONLINEAR_SUNDIALS

#define USE_AUTODIFF
//#define USE_AUTODIFF_ASMJIT
//#define USE_AUTODIFF_EXPRESSION_TEMPLATES

//#define USE_MPI //include mpi for mpi functions
//#define USE_MPI_P2P //use (probably) more effective mpi-2 algorithms
//#define USE_MPI_FILE //use MPI_File_xxx functionality
//#define USE_MPI2 //set of your version produce warnings

//#define USE_INT64
#define USE_FP64
#endif //INMOST_OPTIONS_CMAKE_INCLUDED


// a very expensive check for debug purposes, 
// when you release marker checks all the elements
// that no element is marked by released marker
// in Mesh::Init function change two variables:
// check_shared_mrk  - check shared markers.
// check_private_mrk - check private markers.
//#define CHECKS_MARKERS

// use additional sets to store elements for parallel
// exchanges, otherwise it will recompute those elements
// which is quite expensive
#define USE_PARALLEL_STORAGE

// output xml files for debugging of parallel algorithms
// search for style.xsl within examples for comfortable
// view of generated xml files
//#define USE_PARALLEL_WRITE_TIME

// this will revert Mesh::PrepareReceiveInner to always
// use MPI point to point functionality disregarding problem type
// read comments in Mesh::SetParallelStrategy for more info
// USE_MPI_P2P must be enabled for feature to work
//#define PREFFER_MPI_P2P

// this will write out shared set of elements that lay between processors in GMV format
//#define DEBUG_COMPUTE_SHARED_SKIN_SET

// this controls how sparse data will be allocated,
// when difinition is not commented memory will be always consumed
// by data structure needed to support sparse data,
// otherwise data structure will be allocated only when
// sparse data is presenet
//#define LAZY_SPARSE_ALLOCATION

// this will force array class to use 12 bytes instead of 16 bytes
#define PACK_ARRAY

#define __INLINE inline

#if defined(USE_OMP)
#include <omp.h>
#endif
#if defined(USE_MPI)
#define NOMINMAX
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>
#if !defined(MSMPI_VER) && !defined(MPIO_INCLUDE) && defined(USE_MPI_FILE) && !defined(OMPI_PROVIDE_MPI_FILE_INTERFACE)
//#include <mpio.h> //some versions of MPI doesn't include that
#endif
#endif

#if defined(USE_OPENCL)
#ifdef __APPLE__
#include "OpenCL/opencl.h"
#else
#include "CL/cl.h"
#endif
#endif

#include <string>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

#if defined(USE_OMP)
#define OMP_THREAD omp_get_thread_num()
#define MAX_THREADS omp_get_max_threads()
#else //USE_OMP
#define OMP_THREAD 0
#define MAX_THREADS 1
#endif //USE_OMP


#if defined(min)
#undef min
#endif
#if defined(max)
#undef max
#endif

#if !defined(USE_MPI)
#define INMOST_MPI_Request     int
#define INMOST_MPI_Type        int
#define INMOST_MPI_Comm        int
#define INMOST_MPI_Group       int
#define INMOST_MPI_COMM_WORLD  0
#define INMOST_MPI_COMM_SELF   0
#define INMOST_MPI_BYTE        0
#define INMOST_MPI_INT         0
#define INMOST_MPI_INT64_T     0
#define INMOST_MPI_UINT64_T    0
#define INMOST_MPI_DOUBLE      0
#define INMOST_MPI_FLOAT       0
#define INMOST_MPI_UNSIGNED    0
#define INMOST_MPI_UNSIGNEDL   0
#define INMOST_MPI_UNSIGNEDLL  0
#define INMOST_MPI_Win         int
#define INMOST_MPI_DATATYPE_NULL 0
#define INMOST_MPI_GROUP_EMPTY 0
#else
#define INMOST_MPI_Request     MPI_Request
#define INMOST_MPI_Type        MPI_Datatype
#define INMOST_MPI_Comm        MPI_Comm
#define INMOST_MPI_Group       MPI_Group
#define INMOST_MPI_COMM_WORLD  MPI_COMM_WORLD
#define INMOST_MPI_COMM_SELF   MPI_COMM_SELF
#define INMOST_MPI_BYTE        MPI_BYTE
#define INMOST_MPI_INT         MPI_INT
#define INMOST_MPI_INT64_T     MPI_INT64_T
#define INMOST_MPI_UINT64_T    MPI_UINT64_T
#define INMOST_MPI_DOUBLE      MPI_DOUBLE
#define INMOST_MPI_FLOAT       MPI_FLOAT
#define INMOST_MPI_UNSIGNED    MPI_UNSIGNED
#define INMOST_MPI_UNSIGNEDL   MPI_UNSIGNED_LONG
#define INMOST_MPI_UNSIGNEDLL  MPI_UNSIGNED_LONG_LONG 
#define INMOST_MPI_Win         MPI_Win
#define INMOST_MPI_DATATYPE_NULL MPI_DATATYPE_NULL
#define INMOST_MPI_GROUP_EMPTY MPI_GROUP_EMPTY
#endif

#if !defined(USE_OMP)
#define INMOST_OMP_LOCK_T      int
#else
#define INMOST_OMP_LOCK_T      omp_lock_t
#endif


#define INMOST_MPI_SIZE           int //in case MPI standard changes and compiler gives tons of warnings

#define INMOST_DATA_BULK_TYPE     unsigned char //this should be one byte long
#define INMOST_MPI_DATA_BULK_TYPE     INMOST_MPI_BYTE

#if defined(USE_FP64)
#define INMOST_DATA_REAL_TYPE     double
#define INMOST_MPI_DATA_REAL_TYPE     INMOST_MPI_DOUBLE
#else //USE_FP64
#define INMOST_DATA_REAL_TYPE     float
#define INMOST_MPI_DATA_REAL_TYPE     INMOST_MPI_FLOAT
#endif //USE_FP64


#if defined(USE_INT64)
#define INMOST_DATA_INTEGER_TYPE   int64_t
#define INMOST_DATA_ENUM_TYPE      uint64_t
#define INMOST_DATA_BIG_ENUM_TYPE  uint64_t

#ifndef UINT64_MAX
#define UINT64_MAX ULLONG_MAX
#endif //UINT64_MAX

#define ENUMUNDEF                 UINT64_MAX
#define BIGENUMUNDEF              UINT64_MAX


#define INMOST_MPI_DATA_INTEGER_TYPE   INMOST_MPI_INT64_T
#define INMOST_MPI_DATA_ENUM_TYPE      INMOST_MPI_UINT64_T
#define INMOST_MPI_DATA_BIG_ENUM_TYPE  INMOST_MPI_UINT64_T
#else //USE_INT64
#define INMOST_DATA_INTEGER_TYPE    int           
#define INMOST_DATA_ENUM_TYPE       unsigned int
#define INMOST_DATA_BIG_ENUM_TYPE   unsigned long long

#define ENUMUNDEF                     UINT_MAX
#define BIGENUMUNDEF                  ULLONG_MAX

#define INMOST_MPI_DATA_INTEGER_TYPE   INMOST_MPI_INT
#define INMOST_MPI_DATA_ENUM_TYPE      INMOST_MPI_UNSIGNED
#define INMOST_MPI_DATA_BIG_ENUM_TYPE  INMOST_MPI_UNSIGNEDLL
#endif //USE_INT64

/// Cross-platform timer that return current time in seconds.
/// The timer is similar to MPI_Wtime() and omp_get_wtime() but is independent on both flags USE_MPI and USE_OMP.
///
/// double seconds = Timer();
/// ...some code...
/// seconds = Timer() - seconds;
/// std::cout << "Time spent is " << seconds << " seconds." << std::endl;
double Timer();

namespace INMOST
{
	/// Types of errors may occur in INMOST.
	/// All of these error are fatal ones.
	/// If error is detected then "throw" exception is generated.
	/// The names of the error type are very intuitive and self-explained ones.
	/// Use "try{...} catch{...}" blocks to implement exception handling.
	enum ErrorType 
	{
		/// The list of errors connected to mesh consistency.
		Failure=100,
		NoTagPosition,
		WrongDataType,
		WrongElementType,
		BadTag,
		BadBulkType,
		NoData,
		TagNotInitialized,
		TagNotFound, 
		TagExists, 
		TagForOtherMesh,
		ElementForOtherMesh,
		ImpossibleConn,
		NoElementType,
		NoMultiElement,
		NoMeshElement,
		NoEsetElement,
		DimensionIsFixed,
		NullInElementSet,		
		NoSpaceForMarker,
		ElementBelongsToNobody,

		/// The list of general type errors.
		BadFileName,
		BadFile,		
		CorruptedIerarchy,
		CorruptedOrdering,
		IterForOtherMesh,
		UndefinedBehaviorInGeometry,
		DimensionIsNotSupportedByGeometry,
		NoParallelMode,
		BadParameter,
		TopologyCheckError,
		
		/// The list of errors may occur in the Linear Solver.
		ErrorInSolver = 400,
		DataCorruptedInSolver,
		DifferentCommunicatorInSolver,
		MatrixNotSetInSolver,
		InconsistentSizesInSolver,
		IndexesAreDifferentInSolver,
		PrepareMatrixFirst,
		CannotReusePreconditionerOfDifferentSize,
		SolverNotFound,
		SolverUnsupportedOperation,
		SolverUnknownParameter,
		SolverNonExistentParameters,
		SolverCopyNullException,
		SolverCopyException,
        SolverAssignNullException,
        SolverAssignException,
		
		/// The list of errors may occur in the Partitioner.
		ErrorInPartitioner = 500,
		UnknownWeightSize,
		DistributionTagWasNotFilled,
		
		/// The list of errros that may occur in linear algebra
		MatrixError = 600,
		MatrixCholeskySolveFail,
		MatrixSolveFail,
		MatrixPseudoSolveFail,
		
		/// The very tail of the errors list.
		NotImplemented = 1000,
		Impossible
	};
}

#include "container.hpp"
//#include "io.hpp"


namespace INMOST
{
	template<typename Var>
	class AbstractMatrixReadOnly;

	template<typename Var>
	class AbstractMatrix;

	template<typename Var>
	class MatrixUnit;

	template<typename Var>
	class MatrixDiag;
	
	template<typename Var>
	class SubMatrix;
	
	template<typename Var>
	class ConstSubMatrix;
	
	template<typename Var>
	class BlockOfMatrix;
	
	template<typename Var>
	class ConstBlockOfMatrix;

	template<typename VarA, typename VarB>
	class MatrixSum;

	template<typename VarA, typename VarB>
	class MatrixDifference;

	template<typename Var>
	class MatrixTranspose;

	template<typename VarA, typename VarB>
	class KroneckerProduct;

	template<typename VarA, typename VarB, typename VarR>
	class MatrixMulCoef;

	template<typename VarA, typename VarB, typename VarR>
	class MatrixDivCoef;

	template<typename VarA, typename VarB, typename VarR>
	class MatrixMulShellCoef;

	template<typename VarA, typename VarB, typename VarR>
	class MatrixDivShellCoef;

	template<typename VarA, typename VarB, typename VarR>
	class MatrixMul;
	
	template<typename Var, typename Storage = std::vector<Var> >
	class Matrix;
	
	template<typename Var, typename Storage = std::vector<Var> >
	class SymmetricMatrix;

	class value_reference;
}

#include <cfloat>

__INLINE bool __isnan__(double x) { return x != x; }
__INLINE bool __isinf__(double x) { return fabs(x) > DBL_MAX; }
__INLINE bool __isbad(double x) { return __isnan__(x) || __isinf__(x); }

//in Utils/mem.cpp
size_t getPeakRSS();
size_t getCurrentRSS();

#endif //INMOST_COMMON_INCLUDED
