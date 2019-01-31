cmake_minimum_required (VERSION 2.6)
project (INMOST)


include_directories("Source/Headers")

set(INMOST_MAJOR_VERSION 0)
set(INMOST_MINOR_VERSION 1)
set(INMOST_PATCH_VERSION 0)
set(INMOST_VERSION "${INMOST_MAJOR_VERSION}.${INMOST_MINOR_VERSION}.${INMOST_PATCH_VERSION}")



option(USE_MPI "Compile with MPI support" ON)
option(USE_MPI_P2P "Use MPI point to point functionality, may be faster with hardware support" ON)
option(USE_MPI_FILE "Use MPI extension to work with files, may save a lot of memory" ON)
option(USE_MPI2 "Use MPI-2 extensions, useful if your MPI library warns you to use new functions" ON)
option(USE_OMP "Compile with OpenMP support (experimental)" OFF)
option(USE_OPENCL "Use OpenCL where possible (experimental)" OFF)

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
option(USE_SOLVER_TTSP_FEATURE "Use Trace and Tune Solver Parameters optimization toolkit" OFF)
#option(USE_AUTODIFF_ASMJIT "Use AsmJit for automatic differentiation" OFF)
#option(USE_AUTODIFF_EXPRESSION_TEMPLATES "Use c++ expression templates for automatic differentiation" OFF)

if( MSVC )
option(USE_MT "Use /MT switch instead of /MD switch for libraries generated with Visual Studio" OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
endif( MSVC )

add_subdirectory(Source)
add_definitions(${SOLVER_DEFINITIONS})

add_library(inmost STATIC ${SOURCE} ${HEADER} Source/Solver/ttsp/optimization_parameters.cpp Source/Solver/ttsp/optimization_parameters.h Source/Solver/ttsp/optimizers/bruteforce/ttsp_bruteforce.cpp Source/Solver/ttsp/optimizers/bruteforce/ttsp_bruteforce.h)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_find/")


if( USE_MT )
set(CompilerFlags CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE)
foreach(CompilerFlag ${CompilerFlags}) 
	string(REPLACE "/MD" "/MT" ${CompilerFlag} "${${CompilerFlag}}")
endforeach()
endif( USE_MT )

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
		set(USE_MPI OFF CACHE BOOL "Compile with MPI support" FORCE)
		message("MPI NOT FOUND")
	else()
		include_directories(${MPI_INCLUDE_PATH})
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
		set(USE_OMP OFF CACHE BOOL "Compile with OpenMP support" FORCE)
		message("OpenMP NOT FOUND")
	endif()
endif()

if(USE_PARTITIONER_PARMETIS)
	find_package(ParMETIS)
	if(NOT PARMETIS_FOUND)
		set(USE_PARTITIONER_PARMETIS OFF CACHE BOOL "Use ParMetis partitioner" FORCE)
		message("PARMETIS NOT FOUND")
	else()
		include_directories(${PARMETIS_INCLUDE_DIR})
		include_directories(${METIS_INCLUDE_DIR})
		message("PARMETIS FOUND")
	endif()
endif()

if(USE_SOLVER_METIS)
	find_package(METIS)
	if(NOT METIS_FOUND)
		set(USE_SOLVER_METIS OFF CACHE BOOL "Use METIS for matrix reordering" FORCE)
		message("METIS NOT FOUND")
	else()
		include_directories(${METIS_INCLUDE_DIR})
		message("METIS FOUND")
	endif()
endif()

if(USE_SOLVER_MONDRIAAN)
	find_package(MONDRIAAN)
	if(NOT MONDRIAAN_FOUND)
		set(USE_SOLVER_MONDRIAAN OFF CACHE BOOL "Use Mondriaan for matrix reordering" FORCE)
		message("Mondriaan NOT FOUND")
	else()
		link_directories(${MONDRIAAN_LIBRARY_DIRS})
		include_directories(${MONDRIAAN_INCLUDE_DIRS})
		message("Mondriaan FOUND")
	endif()
endif()


if(USE_PARTITIONER_ZOLTAN)
	find_package(ZOLTAN)
	if(NOT ZOLTAN_FOUND)
		set(USE_PARTITIONER_ZOLTAN OFF CACHE BOOL "Use Zoltan partitioner" FORCE)
		message("ZOLTAN NOT FOUND")
	else()
		#include_directories(${ZOLTAN_INCLUDE_DIR})
		include_directories(${ZOLTAN_INCLUDE_DIRS})
		message("ZOLTAN FOUND")
	endif()
endif()

if(USE_SOLVER_PETSC)
	find_package(PETSc)
	if(NOT PETSC_FOUND)
		set(USE_SOLVER_PETSC OFF CACHE BOOL "Use PETSc solvers" FORCE)
		message("PETSC NOT FOUND")
	else()
		include_directories(${PETSC_INCLUDES})
		message("PETSC FOUND")
		add_definitions(${PETSC_DEFINITIONS})
		#message(${PETSC_LIBRARIES})
	endif()
endif()


if(USE_SOLVER_TRILINOS)
	set(CMAKE_PREFIX_PATH ${TRILINOS_PATH} ${CMAKE_PREFIX_PATH})
	find_package(Trilinos PATHS ${TRILINOS_PATH}/lib/cmake/Trilinos ${TRILINOS_PATH})
	if(NOT Trilinos_FOUND)
		set(USE_SOLVER_TRILINOS OFF CACHE BOOL "Use Trilinos solvers" FORCE)
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


		message("Trilinos FOUND")
	endif()
endif()


if(USE_SOLVER_SUPERLU)
	find_package(SUPERLU)
	if(NOT SUPERLU_FOUND)
		set(USE_SOLVER_SUPERLU OFF CACHE BOOL "Use SuperLU solver" FORCE)
		MESSAGE("SUPERLU NOT FOUND")
	else()
		MESSAGE("\nFound SuperLU!  Here are the details: ")
		MESSAGE("INCLUDES:  ${SUPERLU_INCLUDES}")
		MESSAGE("LIBRARIES: ${SUPERLU_LIBRARIES}")

		include_directories(${SUPERLU_INCLUDES})
	endif()
endif()

if(USE_OPENCL)
	find_package(OpenCL)
	if(OPENCL_FOUND)
		include_directories(${OPENCL_INCLUDE_DIRS})
		link_directories(${OpenCL_LIBRARY})
	else()
		set(USE_OPENCL OFF CACHE BOOL "Use OpenCL where possible (experimental)" FORCE)
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
	message("Adding flag for exception handling in Visual Studio that is sometimes missing for some reason")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
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
                           Source/Headers/inmost_residual.h
                           Source/Headers/inmost_model.h
                           Source/Headers/inmost_operator.h
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
                           Source/Headers/inmost_block_variable.h
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
	"${PROJECT_SOURCE_DIR}/Source/Headers/inmost_residual.h"
	"${PROJECT_SOURCE_DIR}/Source/Headers/inmost_operator.h"
	"${PROJECT_SOURCE_DIR}/Source/Headers/inmost_model.h"
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
	"${PROJECT_SOURCE_DIR}/Source/Headers/inmost_block_variable.h"
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
