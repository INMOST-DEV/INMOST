#pragma once
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
#define USE_AUTODIFF_ASMJIT
//#define USE_AUTODIFF_OPENCL
//#define USE_AUTODIFF_EXPRESSION_TEMPLATES

//#define USE_MPI //include mpi for mpi functions
//#define USE_MPI2 //use (probably) more effective mpi-2 algorithms
#endif //INMOST_OPTIONS_CMAKE_INCLUDED



#define USE_QSORT //use qsort instead of std::sort
#define USE_PARALLEL_STORAGE
#define USE_PARALLEL_WRITE_TIME


#define USE_COMPARE CompareElementsPointer
//#define USE_COMPARE CompareElementsCentroid
//#define USE_COMPARE CompareElementsUnique
//#define USE_COMPARE CompareElementsHybrid


#define __INLINE inline

#if defined(USE_OMP)
#include <omp.h>
#endif
#if defined(USE_MPI)
#include <mpi.h>
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


long double Timer();

namespace INMOST
{
	enum ErrorType 
	{
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
		
		//solver
		ErrorInSolver = 400,
		DataCorruptedInSolver,
		DifferentCommunicatorInSolver,
		MatrixNotSetInSolver,
		InconsistentSizesInSolver,
		IndexesAreDifferentInSolver,
		PrepareMatrixFirst,
		CannotReusePreconditionerOfDifferentSize,
		
		//partitioner
		ErrorInPartitioner = 500,
		UnknownWeightSize,
		DistributionTagWasNotFilled,
		
		NotImplemented = 1000,
		Impossible
	};
}

#include "container.hpp"
#include "io.hpp"


#endif //INMOST_COMMON_INCLUDED
