cmake_minimum_required (VERSION 2.6)
project (INMOST)
set(SOURCE solver.cpp 
           solver_ani.cpp 
           solver_petsc.cpp 
           partitioner.cpp 
           algorithm.cpp 
           geometry.cpp 
           iterator.cpp
           storage.cpp
           eset.cpp
           mesh_file.cpp
           timer.cpp
           face.cpp
           edge.cpp
           mesh.cpp
           node.cpp
           cell.cpp
           tag.cpp
           element.cpp
           mesh_parallel.cpp
           modify.cpp
		   autodiff.cpp)
set(HEADER inmost.h
           inmost_options_cmake.h
           inmost_common.h
           inmost_mesh.h
           inmost_solver.h
           inmost_partitioner.h
           inmost_autodiff.h
           container.hpp
           io.hpp
           solver_ilu2.hpp
           solver_ddpqiluc2.hpp
           solver_bcgsl.hpp
           solver_prototypes.hpp)




add_library(inmost STATIC ${SOURCE} ${HEADER})

option(USE_MPI "Compile with MPI support" ON)
option(USE_OMP "Compile with OpenMP support (experimental)" OFF)

option(USE_MESH "Compile mesh capabilities" ON)
option(USE_SOLVER "Compile solver capabilities" ON)
option(USE_PARTITIONER "Compile partitioner capabilities" ON)
option(USE_AUTODIFF "Compile automatic differentiation capabilities" ON)
option(TEST_FORTRAN_ANI3D "Test for fortran availibility to compile ANI3D lib" OFF)
option(COMPILE_EXAMPLES "Compile examples" OFF)
option(COMPILE_PROJECTS "Compile projects" OFF)
option(COMPILE_TESTS "Compile some tests" OFF)

option(USE_PARTITIONER_PARMETIS "Use ParMetis partitioner" OFF)
option(USE_PARTITIONER_ZOLTAN "Use Zoltan partitioner" OFF)
option(USE_SOLVER_PETSC "Use PETSc solver" OFF)
option(USE_AUTODIFF_OPENCL "Use OpenCL for automatic differentiation (under work)" OFF)
option(USE_AUTODIFF_ASMJIT "Use AsmJit for automatic differentiation" OFF)
option(USE_AUTODIFF_EXPRESSION_TEMPLATES "Use c++ expression templates for automatic differentiation" OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_find/")


if(TEST_FORTRAN_ANI3D)
	add_subdirectory(ILU2)
	if(ANI3D_AVAILIBLE)
		option(USE_SOLVER_ANI "Compile with ANI3D sequential solvers" OFF)
	else()
		option(USE_SOLVER_ANI "Compile with ANI3D sequential solvers" ON)
	endif()
endif()


if(USE_MPI)
	find_package(MPI)
	if(NOT MPI_FOUND)
		option(USE_MPI "Use Message Passing Interface" OFF)
		option(USE_MPI2 "Use MPI-2 functions" OFF)
		message("MPI NOT FOUND")
	else()
		include_directories(${MPI_INCLUDE_PATH})
		option(USE_MPI "Use Message Passing Interface" ON)
		option(USE_MPI2 "Use MPI-2 functions" ON)
		set_target_properties(inmost PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
		message("MPI FOUND")
	endif()
endif()

if(USE_OMP)
	find_package(OpenMP)
	if (OPENMP_FOUND)
		set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	else()
		option(USE_OMP "Compile with OpenMP support (experimental)" OFF)
	endif()
endif()

if(USE_PARTITIONER_PARMETIS)
	find_package(ParMETIS)
	if(NOT PARMETIS_FOUND)
		option(USE_PARTITIONER_PARMETIS "Use ParMetis partitioner" OFF)
		message("PARMETIS NOT FOUND")
	else()
		include_directories(${PARMETIS_INCLUDE_DIR})
		option(USE_PARTITIONER_PARMETIS "Use ParMetis partitioner" ON)
		message("PARMETIS FOUND")
	endif()
endif()


if(USE_PARTITIONER_ZOLTAN)
	find_package(ZOLTAN)
	if(NOT ZOLTAN_FOUND)
		option(USE_PARTITIONER_ZOLTAN "Use Zoltan partitioner" OFF)
		message("ZOLTAN NOT FOUND")
	else()
		include_directories(${ZOLTAN_INCLUDE_DIRS})
		option(USE_PARTITIONER_ZOLTAN "Use Zoltan partitioner" ON)
		message("ZOLTAN FOUND")
	endif()
endif()

if(USE_SOLVER_PETSC)
	find_package(PETSc)
	if(NOT PETSC_FOUND)
		option(USE_SOLVER_PETSC "Use PETSc solver" OFF)
		message("PETSC NOT FOUND")
	else()
		include_directories(${PETSC_INCLUDES})
		option(USE_SOLVER_PETSC "Use PETSc solver" ON)
		message("PETSC FOUND")
		add_definitions(${PETSC_DEFINITIONS})
		#message(${PETSC_LIBRARIES})
	endif()
endif()


if(USE_AUTODIFF_OPENCL)
	find_package(OpenCL)
	if(OPENCL_FOUND)
		option(USE_AUTODIFF_OPENCL "Use OpenCL for automatic differentiation" ON)
		include_directories(${OPENCL_INCLUDE_DIRS})
	else()
		option(USE_AUTODIFF_OPENCL "Use OpenCL for automatic differentiation" OFF)
		message("OpenCL not found")
	endif()
endif()

if(USE_AUTODIFF_ASMJIT)
	Set(ASMJIT_STATIC TRUE)
	include("AsmJit/CMakeLists.txt")	
	include_directories("${ASMJIT_INC_DIR}/asmjit")
	target_link_libraries(inmost asmjit)
endif()

if(USE_AUTODIFF_EXPRESSION_TEMPLATES)
	if(MSVC)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
	endif()
endif()


configure_file("inmost_options_cmake.h" "${PROJECT_BINARY_DIR}/inmost_options.h")
include_directories("${PROJECT_BINARY_DIR}")

if(COMPILE_EXAMPLES)
add_subdirectory(examples)
endif(COMPILE_EXAMPLES)

if(COMPILE_PROJECTS)
add_subdirectory(projects)
endif(COMPILE_PROJECTS)

if(COMPILE_TESTS)
add_subdirectory(tests)
endif(COMPILE_TESTS)

