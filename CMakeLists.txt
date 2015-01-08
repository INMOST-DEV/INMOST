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
           earray.cpp
           comparator.cpp
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
option(USE_MPI_P2P "Use MPI point to point functionality, may be faster with hardware support" ON)
option(USE_MPI_FILE "Use MPI extension to work with files, may save a lot of memory" ON)
option(USE_MPI2 "Use MPI-2 extensions, useful if your MPI library warns you to use new functions" ON)
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
		set(USE_MPI OFF)
		message("MPI NOT FOUND")
	else()
		include_directories(${MPI_INCLUDE_PATH})
		set(USE_MPI ON)
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
		set(USE_OMP OFF)
	endif()
endif()

if(USE_PARTITIONER_PARMETIS)
	find_package(ParMETIS)
	if(NOT PARMETIS_FOUND)
		set(USE_PARTITIONER_PARMETIS OFF)
		message("PARMETIS NOT FOUND")
	else()
		include_directories(${PARMETIS_INCLUDE_DIR})
		include_directories(${METIS_INCLUDE_DIR})
		set(USE_PARTITIONER_PARMETIS ON)
		message("PARMETIS FOUND")
	endif()
endif()


if(USE_PARTITIONER_ZOLTAN)
	find_package(ZOLTAN)
	if(NOT ZOLTAN_FOUND)
		set(USE_PARTITIONER_ZOLTAN OFF)
		message("ZOLTAN NOT FOUND")
	else()
		include_directories(${ZOLTAN_INCLUDE_DIR})
		include_directories(${ZOLTAN_INCLUDE_DIRS})
		set(USE_PARTITIONER_ZOLTAN ON)
		message("ZOLTAN FOUND")
	endif()
endif()

if(USE_SOLVER_PETSC)
	find_package(PETSc)
	if(NOT PETSC_FOUND)
		set(USE_SOLVER_PETSC OFF)
		message("PETSC NOT FOUND")
	else()
		include_directories(${PETSC_INCLUDES})
		set(USE_SOLVER_PETSC ON)
		message("PETSC FOUND")
		add_definitions(${PETSC_DEFINITIONS})
		#message(${PETSC_LIBRARIES})
	endif()
endif()


if(USE_AUTODIFF_OPENCL)
	find_package(OpenCL)
	if(OPENCL_FOUND)
		set(USE_AUTODIFF_OPENCL ON)
		include_directories(${OPENCL_INCLUDE_DIRS})
	else()
		set(USE_AUTODIFF_OPENCL OFF)
		message("OpenCL not found")
	endif()
endif()

if(USE_AUTODIFF_ASMJIT)
	set(ASMJIT_STATIC TRUE)
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

