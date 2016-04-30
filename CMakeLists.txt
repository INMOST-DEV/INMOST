cmake_minimum_required (VERSION 2.6)
project (INMOST)

add_subdirectory(Source)
add_definitions(${SOLVER_DEFINITIONS})
include_directories("Source/Headers")

set(INMOST_MAJOR_VERSION 0)
set(INMOST_MINOR_VERSION 1)
set(INMOST_PATCH_VERSION 0)
set(INMOST_VERSION "${INMOST_MAJOR_VERSION}.${INMOST_MINOR_VERSION}.${INMOST_PATCH_VERSION}")


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
option(USE_NONLINEAR "Compile nonlinear solver capabilities" ON)
option(TEST_FORTRAN_ANI3D "Test for fortran availibility to compile ANI3D lib" OFF)
option(COMPILE_EXAMPLES "Compile examples" OFF)
option(COMPILE_TESTS "Compile some tests" OFF)

option(USE_PARTITIONER_PARMETIS "Use ParMetis partitioner" OFF)
option(USE_PARTITIONER_ZOLTAN "Use Zoltan partitioner" OFF)
option(USE_SOLVER_METIS "Use METIS for matrix reordering" OFF)
option(USE_SOLVER_MONDRIAAN "Use Mondriaan for matrix reordering" OFF)
option(USE_SOLVER_PETSC "Use PETSc solvers" OFF)
option(USE_SOLVER_TRILINOS "Use Trilinos solvers" OFF)
option(USE_SOLVER_SUPERLU "Use SuperLU solver" OFF)
#option(USE_AUTODIFF_OPENCL "Use OpenCL for automatic differentiation (under work)" OFF)
#option(USE_AUTODIFF_ASMJIT "Use AsmJit for automatic differentiation" OFF)
#option(USE_AUTODIFF_EXPRESSION_TEMPLATES "Use c++ expression templates for automatic differentiation" OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_find/")


if(TEST_FORTRAN_ANI3D)
	include(CheckLanguage)
	check_language(Fortran)
	find_package(BLAS)
	if(NOT CMAKE_Fortran_COMPILER OR NOT BLAS_FOUND)
		message(STATUS "No Fortran or BLAS support")
		option(USE_SOLVER_ANI "Compile with ANI3D sequential solvers" OFF)
	else()
		enable_language(Fortran)
		add_subdirectory(ILU2)
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
		message("OpenMP FOUND")
	else()
		set(USE_OMP OFF)
		message("OpenMP NOT FOUND")
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

if(USE_SOLVER_METIS)
	find_package(METIS)
	if(NOT METIS_FOUND)
		set(USE_SOLVER_METIS OFF)
		message("METIS NOT FOUND")
	else()
		include_directories(${METIS_INCLUDE_DIR})
		set(USE_SOLVER_METIS ON)
		message("METIS FOUND")
	endif()
endif()

if(USE_SOLVER_MONDRIAAN)
	find_package(MONDRIAAN)
	if(NOT MONDRIAAN_FOUND)
		set(USE_SOLVER_MONDRIAAN OFF)
		message("Mondriaan NOT FOUND")
	else()
		link_directories(${MONDRIAAN_LIBRARY_DIRS})
		include_directories(${MONDRIAAN_INCLUDE_DIRS})
		set(USE_SOLVER_MONDRIAAN ON)
		message("Mondriaan FOUND")
	endif()
endif()


if(USE_PARTITIONER_ZOLTAN)
	find_package(ZOLTAN)
	if(NOT ZOLTAN_FOUND)
		set(USE_PARTITIONER_ZOLTAN OFF)
		message("ZOLTAN NOT FOUND")
	else()
		#include_directories(${ZOLTAN_INCLUDE_DIR})
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
		include_directories(${PETSC_INCLUDE_DIR})
		set(USE_SOLVER_PETSC ON)
		message("PETSC FOUND")
		add_definitions(${PETSC_DEFINITIONS})
		#message(${PETSC_LIBRARIES})
	endif()
endif()


if(USE_SOLVER_TRILINOS)
	set(CMAKE_PREFIX_PATH ${TRILINOS_PATH} ${CMAKE_PREFIX_PATH})
	find_package(Trilinos PATHS ${TRILINOS_PATH}/lib/cmake/Trilinos ${TRILINOS_PATH})
	if(NOT Trilinos_FOUND)
		set(USE_SOLVER_TRILINOS OFF)
		message("Trilinos NOT FOUND")
	else()

		MESSAGE("\nFound Trilinos!  Here are the details: ")
		MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
		MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
		MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
		MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
		MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
		MESSAGE("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
		MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
		MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
		MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
		MESSAGE("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
		MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
		MESSAGE("End of Trilinos details\n")

		include_directories(${Trilinos_INCLUDE_DIRS})
		include_directories(${Trilinos_TPL_INCLUDE_DIRS})
		link_directories(${Trilinos_LIBRARY_DIRS})
		link_directories(${Trilinos_TPL_LIBRARY_DIRS})


		set(USE_SOLVER_TRILINOS ON)
		message("Trilinos FOUND")
	endif()
endif()


if(USE_SOLVER_SUPERLU)
	find_package(SUPERLU)
	if(NOT SUPERLU_FOUND)
		set(USE_SOLVER_SUPERLU OFF)
		MESSAGE("SUPERLU NOT FOUND")
	else()
		MESSAGE("\nFound SuperLU!  Here are the details: ")
		MESSAGE("INCLUDES:  ${SUPERLU_INCLUDES}")
		MESSAGE("LIBRARIES: ${SUPERLU_LIBRARIES}")

		include_directories(${SUPERLU_INCLUDES})
		set(USE_SOLVER_SUPERLU ON)
		
	endif()
endif()

if(USE_AUTODIFF_OPENCL)
	find_package(OpenCL)
	if(OPENCL_FOUND)
		set(USE_AUTODIFF_OPENCL ON)
		include_directories(${OPENCL_INCLUDE_DIRS})
		link_directories(${PETSC_LIBRARY_DIRS})
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


if(MSVC)
  if(USE_SOLVER_TRILINOS)
    message("Putting workaround for Visual Studio that allow to use Trilinos Release libraries in Debug mode")
    message("Note that this workaround may affect your debugging experience, you may want to debug without Trilinos")
    add_definitions(-D_ITERATOR_DEBUG_LEVEL=0)
  endif()
endif()

configure_file("Source/Headers/inmost_options_cmake.h" "${PROJECT_BINARY_DIR}/inmost_options.h")
include_directories("${PROJECT_BINARY_DIR}")

target_include_directories(inmost PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_include_directories(inmost PUBLIC $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>)
target_include_directories(inmost PUBLIC $<INSTALL_INTERFACE:include>)

if(COMPILE_EXAMPLES)
add_subdirectory(Examples)
endif(COMPILE_EXAMPLES)


if(COMPILE_TESTS)
enable_testing()
include(CTest)
add_subdirectory(Tests)
endif(COMPILE_TESTS)


set(INMOST_INSTALL_HEADERS Source/Headers/inmost.h
                           Source/Headers/inmost_autodiff.h
                           Source/Headers/inmost_common.h
                           Source/Headers/inmost_data.h
                           Source/Headers/inmost_dense.h
                           Source/Headers/inmost_expression.h
                           Source/Headers/inmost_mesh.h
                           Source/Headers/inmost_nonlinear.h
                           Source/Headers/inmost_partitioner.h
                           Source/Headers/inmost_solver.h
                           Source/Headers/inmost_sparse.h
                           Source/Headers/inmost_xml.h
                           Source/Headers/inmost_variable.h
                           Source/Headers/container.hpp)



include(CPack)
export(TARGETS inmost FILE inmost-targets.cmake)
export(PACKAGE inmost)
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/Source/Headers" "${PROJECT_BINARY_DIR}")
configure_file(inmost-config.cmake.in "${PROJECT_BINARY_DIR}/inmost-config.cmake" @ONLY)
set(CONF_INCLUDE_DIRS "\${inmost_DIR}/include")
configure_file(inmost-config.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/inmost-config.cmake" @ONLY)
configure_file(inmost-config-version.cmake.in "${PROJECT_BINARY_DIR}/inmost-config-version.cmake" @ONLY)



set_property(TARGET inmost PROPERTY PUBLIC_HEADER
	"${PROJECT_BINARY_DIR}/inmost_options.h"
	"${PROJECT_SOURCE_DIR}/Source/Headers/inmost.h"
	"${PROJECT_SOURCE_DIR}/Source/Headers/inmost_autodiff.h"
	"${PROJECT_SOURCE_DIR}/Source/Headers/inmost_common.h"
	"${PROJECT_SOURCE_DIR}/Source/Headers/inmost_data.h"
	"${PROJECT_SOURCE_DIR}/Source/Headers/inmost_dense.h"
	"${PROJECT_SOURCE_DIR}/Source/Headers/inmost_expression.h"
	"${PROJECT_SOURCE_DIR}/Source/Headers/inmost_mesh.h"
	"${PROJECT_SOURCE_DIR}/Source/Headers/inmost_nonlinear.h"
	"${PROJECT_SOURCE_DIR}/Source/Headers/inmost_partitioner.h"
	"${PROJECT_SOURCE_DIR}/Source/Headers/inmost_solver.h"
	"${PROJECT_SOURCE_DIR}/Source/Headers/inmost_sparse.h"
	"${PROJECT_SOURCE_DIR}/Source/Headers/inmost_variable.h"
	"${PROJECT_SOURCE_DIR}/Source/Headers/inmost_xml.h"
	"${PROJECT_SOURCE_DIR}/Source/Headers/container.hpp")

install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/inmost-config.cmake"
  "${PROJECT_BINARY_DIR}/inmost-config-version.cmake"
  DESTINATION "${CMAKE_INSTALL_PREFIX}")

install(TARGETS inmost EXPORT inmost-targets
        ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
	PUBLIC_HEADER DESTINATION include)


install(EXPORT inmost-targets DESTINATION cmake)
