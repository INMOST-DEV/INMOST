
#ifndef INMOST_COMMON_INCLUDED
#define INMOST_COMMON_INCLUDED


#include "inmost_options.h"
#if !defined(INMOST_OPTIONS_CMAKE_INCLUDED)
//no options from cmake -- define minimum
//#define USE_OMP

#define USE_MESH

#define USE_PARTITIONER
//#define USE_PARTITIONER_ZOLTAN
//#define USE_PARTITIONER_PARMETIS

#define USE_SOLVER
//#define USE_SOLVER_PETSC
//#define USE_SOLVER_ANI

#define USE_AUTODIFF
//#define USE_AUTODIFF_ASMJIT
//#define USE_AUTODIFF_OPENCL
//#define USE_AUTODIFF_EXPRESSION_TEMPLATES

//#define USE_MPI //include mpi for mpi functions
//#define USE_MPI_P2P //use (probably) more effective mpi-2 algorithms
//#define USE_MPI_FILE //use MPI_File_xxx functionality
//#define USE_MPI2 //set of your version produce warnings
#endif //INMOST_OPTIONS_CMAKE_INCLUDED


// a very expensive check for debug purposes, 
// when you release marker checks all the elements
// that no element is marked by released marker
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
#include <mpi.h>
#if !defined(MSMPI_VER) && !defined(MPIO_INCLUDE)
#include <mpio.h> //some versions of MPI doesn't include that
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
#include <omp.h>
#endif


#if !defined(USE_MPI)
#define INMOST_MPI_Request     int
#define INMOST_MPI_Type        int
#define INMOST_MPI_Comm        int
#define INMOST_MPI_COMM_WORLD  0
#define INMOST_MPI_BYTE        0
#define INMOST_MPI_INT         0
#define INMOST_MPI_DOUBLE      0
#define INMOST_MPI_UNSIGNED    0
#define INMOST_MPI_Win         int
#else
#define INMOST_MPI_Request     MPI_Request
#define INMOST_MPI_Type        MPI_Datatype
#define INMOST_MPI_Comm        MPI_Comm
#define INMOST_MPI_COMM_WORLD  MPI_COMM_WORLD
#define INMOST_MPI_BYTE        MPI_BYTE
#define INMOST_MPI_INT         MPI_INT
#define INMOST_MPI_DOUBLE      MPI_DOUBLE
#define INMOST_MPI_UNSIGNED    MPI_UNSIGNED
#define INMOST_MPI_Win         MPI_Win
#endif

#define INMOST_MPI_SIZE           int //in case MPI standard changes and compiler gives tons of warning

#define INMOST_DATA_INTEGER_TYPE  int           
#define INMOST_DATA_REAL_TYPE     double        
#define INMOST_DATA_BULK_TYPE     unsigned char //this should be one byte long

#define INMOST_MPI_DATA_INTEGER_TYPE  INMOST_MPI_INT
#define INMOST_MPI_DATA_REAL_TYPE     INMOST_MPI_DOUBLE
#define INMOST_MPI_DATA_BULK_TYPE     INMOST_MPI_BYTE



#define INMOST_DATA_ENUM_TYPE       unsigned int
#define ENUMUNDEF                 UINT_MAX
#define INMOST_DATA_BIG_ENUM_TYPE   unsigned int
#define BIGENUMUNDEF              UINT_MAX

#define INMOST_MPI_DATA_ENUM_TYPE      INMOST_MPI_UNSIGNED
#define INMOST_MPI_DATA_BIG_ENUM_TYPE  INMOST_MPI_UNSIGNED


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
		
		/// The list of errors may occur in the Partitioner.
		ErrorInPartitioner = 500,
		UnknownWeightSize,
		DistributionTagWasNotFilled,
		
		/// The very tail of the errors list.
		NotImplemented = 1000,
		Impossible
	};
}

#include "container.hpp"
#include "io.hpp"


#endif //INMOST_COMMON_INCLUDED
