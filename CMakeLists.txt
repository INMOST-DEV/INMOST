cmake_minimum_required (VERSION 2.6)
project (INMOST)


include_directories("Source/Headers")

set(INMOST_MAJOR_VERSION 0)
set(INMOST_MINOR_VERSION 3)
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
option(USE_PARTITIONER_PARMETIS_DOWNLOAD "Attempts to download ParMetis if library not found" OFF)
option(USE_PARTITIONER_METIS "Use METIS partitioner" OFF)
option(USE_PARTITIONER_ZOLTAN "Use Zoltan partitioner" OFF)
option(USE_SOLVER_METIS "Use METIS for matrix reordering" OFF)
option(USE_SOLVER_MONDRIAAN "Use Mondriaan for matrix reordering" OFF)
option(USE_SOLVER_PETSC "Use PETSc solvers" OFF)
option(USE_SOLVER_TRILINOS "Use Trilinos solvers" OFF)
option(USE_SOLVER_SUPERLU "Use SuperLU solver" OFF)
option(USE_SOLVER_SUPERLU_DOWNLOAD "Attempts to download SuperLU solver if library not found" OFF)
option(USE_OPTIMIZER "Use Optimization toolkit" OFF)
option(USE_OPTIMIZER_BAYESIAN "Compile Optimization module against Limbo library to provide Bayesian optimization algorithm" OFF)
#option(USE_AUTODIFF_ASMJIT "Use AsmJit for automatic differentiation" OFF)
#option(USE_AUTODIFF_EXPRESSION_TEMPLATES "Use c++ expression templates for automatic differentiation" OFF)

set(LIB_DOWNLOAD_PATH "${CMAKE_BINARY_DIR}/Libraries/download")

if( MSVC )
option(USE_MT "Use /MT switch instead of /MD switch for libraries generated with Visual Studio" OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
endif( MSVC )

add_subdirectory(Source)
add_definitions(${SOLVER_DEFINITIONS})

add_library(inmost STATIC ${SOURCE} ${HEADER})

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

if(USE_OPTIMIZER_BAYESIAN)
	include_directories(/Users/bvdmitri/Projects/CXX/limbo/src)
	target_link_libraries(inmost /Users/bvdmitri/Projects/CXX/limbo/build/src/liblimbo.a)

	find_package(Eigen3 REQUIRED)
	include_directories(${EIGEN3_INCLUDE_DIR})

	find_package(NLOPT)
	include_directories(${NLOPT_INCLUDE_DIRS})
	target_link_libraries(inmost ${NLOPT_LIBRARIES})

	find_package(boost REQUIRED)
	include_directories(${Boost_INCLUDE_DIRS})
	target_link_libraries(inmost ${Boost_LIBRARIES})
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
		#set(USE_PARTITIONER_PARMETIS OFF CACHE BOOL "Use ParMetis partitioner" FORCE)
		message("PARMETIS NOT FOUND")

		if( USE_PARTITIONER_PARMETIS_DOWNLOAD )
			include(ExternalProject)
			message("-- Download ParMetis")
			string (REPLACE ";" "|" PASS_MPI "${MPI_INCLUDE_PATH}")
			message("-- Passing MPI include to ParMETIS: '${PASS_MPI}'")
			
			ExternalProject_Add(ParMETIS
								#This one is not good
								#GIT_REPOSITORY "https://github.com/ibaned/parmetis.git"
								GIT_REPOSITORY "https://github.com/scibuilder/parmetis.git"
								GIT_TAG "master"
								UPDATE_DISCONNECTED 1 
								#${EXTERNAL_NOUPDATE}
								PREFIX "${LIB_DOWNLOAD_PATH}"
								#this command is responsible for update of repository from git
								#UPDATE_COMMAND ""
								LIST_SEPARATOR | # Use the alternate list separator
								#for the one from ibaned
								#CMAKE_ARGS "-DMPI_INCLUDE_PATH:PATH='${PASS_MPI}'"
								CMAKE_ARGS  "-DDISABLE_PARMETIS_PROGRAMS:BOOL=TRUE "
											"-DCMAKE_INSTALL_PREFIX=${LIB_DOWNLOAD_PATH} "
											"-DCMAKE_BUILD_TYPE=Release "
											"-DMPI_C_COMPILER=${MPI_C_COMPILER} "
											"-DMPI_CXX_COMPILER=${MPI_CXX_COMPILER} "
											"-DCMAKE_POSITION_INDEPENDENT_CODE=ON "
											"-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} "
											"-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} "
								)
			add_dependencies(inmost ParMETIS)
			set(PARMETIS_INCLUDE_DIR "${LIB_DOWNLOAD_PATH}/include" CACHE PATH "PARMETIS_INCLUDE_DIR is set" FORCE)
			set(METIS_INCLUDE_DIR    "${LIB_DOWNLOAD_PATH}/include" CACHE PATH "METIS_INCLUDE_DIR is set" FORCE)
			if( WIN32 )
				set(PARMETIS_LIBRARIES "${LIB_DOWNLOAD_PATH}/lib/parmetis.lib;${LIB_DOWNLOAD_PATH}/lib/metis.lib")
				set(METIS_LIBRARIES    "${LIB_DOWNLOAD_PATH}/lib/metis.lib")
			else( WIN32 )
				set(PARMETIS_LIBRARIES "${LIB_DOWNLOAD_PATH}/lib/libparmetis.a;${LIB_DOWNLOAD_PATH}/lib/libmetis.a" CACHE PATH "PARMETIS_LIBRARIES is set" FORCE)
				set(METIS_LIBRARIES    "${LIB_DOWNLOAD_PATH}/lib/libmetis.a" CACHE PATH "METIS_LIBRARIES is set" FORCE)
			endif( WIN32 )
			message("-- PARMETIS_INCLUDE_DIR: ${PARMETIS_INCLUDE_DIR}")
			message("-- PARMETIS_LIBRARIES:   ${PARMETIS_LIBRARIES}")
			message("-- METIS_INCLUDE_DIR:    ${METIS_INCLUDE_DIR}")
			message("-- METIS_LIBRARIES:      ${METIS_LIBRARIES}")
			include_directories(${PARMETIS_INCLUDE_DIR})
			include_directories(${METIS_INCLUDE_DIR})
			set(USE_PARTITIONER_METIS ON CACHE BOOL "Use METIS partitioner" FORCE)
		else(USE_PARTITIONER_PARMETIS_DOWNLOAD)
			set(USE_PARTITIONER_PARMETIS OFF CACHE BOOL "Use ParMetis partitioner" FORCE)
		endif(USE_PARTITIONER_PARMETIS_DOWNLOAD)
	else()
		include_directories(${PARMETIS_INCLUDE_DIR})
		include_directories(${METIS_INCLUDE_DIR})
		message("PARMETIS FOUND")
	endif()
endif()

if(USE_PARTITIONER_METIS)
	find_package(METIS)
	if(NOT METIS_FOUND)
		set(USE_SOLVER_METIS OFF CACHE BOOL "Use METIS partitioner" FORCE)
		message("METIS NOT FOUND")
	else()
		include_directories(${METIS_INCLUDE_DIR})
		message("METIS FOUND")
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
	find_package(LAPACK)
	if( NOT LAPACK_FOUND )
		message("-- Download LAPACK")
		ExternalProject_Add(LAPACK
							URL "https://www.netlib.org/clapack/clapack-3.2.1-CMAKE.tgz"
							UPDATE_DISCONNECTED 1
							PREFIX "${LIB_DOWNLOAD_PATH}"
							INSTALL_COMMAND true
							CMAKE_ARGS #"-DCMAKE_INSTALL_PREFIX=${LIB_DOWNLOAD_PATH} " #no install script!
									   "-DMPI_C_COMPILER=${MPI_C_COMPILER} "
									   "-DMPI_CXX_COMPILER=${MPI_CXX_COMPILER} "
									   "-DCMAKE_BUILD_TYPE=Release ")
		if( WIN32 )
			set(LAPACK_LIBRARIES "${LIB_DOWNLOAD_PATH}/src/LAPACK-build/SRC/lapack.lib" CACHE PATH "LAPACK_LIBRARIES is set" FORCE)
			set(BLAS_LIBRARIES "${LIB_DOWNLOAD_PATH}/src/LAPACK-build/BLAS/SRC/blas.lib;${LIB_DOWNLOAD_PATH}/src/LAPACK-build/F2CLIBS/libf2c/f2c.lib" CACHE PATH "LAPACK_LIBRARIES is set" FORCE)
		else( WIN32 )
			set(LAPACK_LIBRARIES "${LIB_DOWNLOAD_PATH}/src/LAPACK-build/SRC/liblapack.a" CACHE PATH "LAPACK_LIBRARIES is set" FORCE)
			set(BLAS_LIBRARIES "${LIB_DOWNLOAD_PATH}/src/LAPACK-build/BLAS/SRC/libblas.a;${LIB_DOWNLOAD_PATH}/src/LAPACK-build/F2CLIBS/libf2c/libf2c.a" CACHE PATH "LAPACK_LIBRARIES is set" FORCE)
		endif( WIN32 )
	endif()
	string (REPLACE ";" "|" PASS_LAPACK_LIBRARIES "${LAPACK_LIBRARIES}")
	string (REPLACE ";" "|" PASS_BLAS_LIBRARIES "${BLAS_LIBRARIES}")
	
	if( NOT USE_MPI )
		find_package(SUPERLU)
		if(NOT SUPERLU_FOUND)
			#set(USE_SOLVER_SUPERLU OFF CACHE BOOL "Use SuperLU solver" FORCE)
			MESSAGE("SUPERLU NOT FOUND")
			if( USE_SOLVER_SUPERLU_DOWNLOAD )
				
				include(ExternalProject)
				message("-- Download SUPERLU")
				#string (REPLACE ";" "|" PASS_MPI "${MPI_INCLUDE_PATH}")
				#message("-- Passing MPI include to ParMETIS: '${PASS_MPI}'")
				
				ExternalProject_Add(SUPERLU
									GIT_REPOSITORY "https://github.com/xiaoyeli/superlu.git"
									GIT_TAG "master"
									UPDATE_DISCONNECTED 1 
									#${EXTERNAL_NOUPDATE}
									PREFIX "${LIB_DOWNLOAD_PATH}"
									#this command is responsible for update of repository from git
									#UPDATE_COMMAND ""
									LIST_SEPARATOR | # Use the alternate list separator
									#for the one from ibaned
									#CMAKE_ARGS "-DMPI_INCLUDE_PATH:PATH='${PASS_MPI}'"
									CMAKE_ARGS  		"-DCMAKE_INSTALL_PREFIX=${LIB_DOWNLOAD_PATH} "
												"-DCMAKE_INSTALL_INCLUDEDIR=${LIB_DOWNLOAD_PATH}/include/superlu "
												"-DBLAS_LIBRARIES=${PASS_BLAS_LIBRARIES} " #should come from find_pacakge(LAPACK)
												"-DCMAKE_BUILD_TYPE=Release "
												"-DMPI_C_COMPILER=${MPI_C_COMPILER} "
												"-DMPI_CXX_COMPILER=${MPI_CXX_COMPILER} "
												"-DCMAKE_POSITION_INDEPENDENT_CODE=ON "
												"-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} "
												"-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} "
									)
				#get internal blas
				#find_package(BLAS)
				#if (NOT BLAS_FOUND)
				#	if( WIN32 )
				#		set(BLAS_LIBRARIES "${LIB_DOWNLOAD_PATH}/src/SUPERLU-build/CBLAS/blas.lib" CACHE FILEPATH "BLAS_LIBRARIES" FORCE)
				#	else( WIN32 )
				#		set(BLAS_LIBRARIES "${LIB_DOWNLOAD_PATH}/src/SUPERLU-build/CBLAS/libblas.a" CACHE FILEPATH "BLAS_LIBRARIES" FORCE)
				#	endif( WIN32 )
				#endif()
				
				add_dependencies(inmost SUPERLU)

				set(SUPERLU_INCLUDES "${LIB_DOWNLOAD_PATH}/include" CACHE PATH "SUPERLU_INCLUDE_DIR is set" FORCE)
				if( WIN32 )
					set(SUPERLU_LIBRARIES "${LIB_DOWNLOAD_PATH}/lib/superlu.lib;${BLAS_LIBRARIES}" CACHE PATH "SUPERLU_LIBRARIES is set" FORCE)
				else( WIN32 )
					include(GNUInstallDirs)
					message("gnu install library dirs ${CMAKE_INSTALL_LIBDIR}")
					set(SUPERLU_LIBRARIES "${LIB_DOWNLOAD_PATH}/${CMAKE_INSTALL_LIBDIR}/libsuperlu.a;${BLAS_LIBRARIES}" CACHE PATH "SUPERLU_LIBRARIES is set" FORCE)
				endif( WIN32 )
				message("-- SUPERLU_INCLUDE_DIR: ${SUPERLU_INCLUDES}")
				message("-- SUPERLU_LIBRARIES:   ${SUPERLU_LIBRARIES}")
			else(USE_SOLVER_SUPERLU_DOWNLOAD)
				set(USE_SOLVER_SUPERLU OFF CACHE BOOL "Use SuperLU solver" FORCE)
			endif(USE_SOLVER_SUPERLU_DOWNLOAD)
		else()
			MESSAGE("\nFound SuperLU!  Here are the details: ")
			MESSAGE("INCLUDES:  ${SUPERLU_INCLUDES}")
			MESSAGE("LIBRARIES: ${SUPERLU_LIBRARIES}")
		endif()
	else()
		find_package(SUPERLU_DIST)
		if(NOT SUPERLU_DIST_FOUND)
			#set(USE_SOLVER_SUPERLU OFF CACHE BOOL "Use SuperLU solver" FORCE)
			MESSAGE("SUPERLU_DIST NOT FOUND")
			
			if( USE_SOLVER_SUPERLU_DOWNLOAD )
				if( USE_MPI AND USE_PARTITIONER_PARMETIS )
					message("-- Download SuperLU_DIST")
					set(C_STANDARD_FLAG "-std=c99")
					string (REPLACE ";" "|" PASS_PARMETIS_LIBRARIES "${PARMETIS_LIBRARIES}")
					message("-- Pass parmetis libraries ${PASS_PARMETIS_LIBRARIES} ")
					message("-- Pass lapack libraries ${PASS_LAPACK_LIBRARIES} ")
					message("-- Pass blas libraries ${PASS_BLAS_LIBRARIES} ")
					ExternalProject_Add(SuperLU_DIST
										GIT_REPOSITORY "https://github.com/xiaoyeli/superlu_dist.git"
										GIT_TAG "master"
										UPDATE_DISCONNECTED 1
										PREFIX "${LIB_DOWNLOAD_PATH}"
										#this command is responsible for update of repository from git
										#UPDATE_COMMAND ""
										LIST_SEPARATOR | # Use the alternate list separator
										CMAKE_ARGS  "-DTPL_PARMETIS_INCLUDE_DIRS=${PARMETIS_INCLUDE_DIR} "
													"-DTPL_PARMETIS_LIBRARIES=${PASS_PARMETIS_LIBRARIES} "
													"-DTPL_LAPACK_LIBRARIES=${PASS_LAPACK_LIBRARIES} "
													"-DBLAS_LIBRARIES=${PASS_BLAS_LIBRARIES} "
													"-Denable_examples=OFF "
													"-Denable_tests=OFF "
													"-Denable_doc=OFF "
													"-DCMAKE_INSTALL_PREFIX=${LIB_DOWNLOAD_PATH} "
													"-DUSE_XSDK_DEFAULTS=OFF "
													"-DUSE_XSDK_DEFAULTS_DEFAULT=OFF "
													"-DXSDK_ENABLE_Fortran=OFF "
													"-DNOFORTRAN=ON "
													"-DCMAKE_C_FLAGS='-DMETIS_EXPORT=' "
													"-DCMAKE_CXX_FLAGS='-DMETIS_EXPORT=' "
													"-DCMAKE_BUILD_TYPE=Release "
													"-DMPI_C_COMPILER=${MPI_C_COMPILER} "
													"-DMPI_CXX_COMPILER=${MPI_CXX_COMPILER} "
													"-DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=ON "
													"-DCMAKE_POSITION_INDEPENDENT_CODE=ON "
										)

					message("-- Linking SuperLU_dist and ParMETIS")
					if( NOT LAPACK_FOUND )
						add_dependencies(SuperLU_DIST LAPACK)
					endif()
					add_dependencies(SuperLU_DIST ParMETIS)
					add_dependencies(inmost SuperLU_DIST)
					add_definitions(-DUSE_SOLVER_SUPERLU_DIST)
					#same as above
					#set(SUPERLU_INCLUDES "${LIB_DOWNLOAD_PATH}/include" CACHE PATH "SUPERLU_INCLUDE_DIR is set" FORCE)
					#extended
					#TODO: properly install lapack
					if( WIN32 )
						set(SUPERLU_LIBRARIES "${SUPERLU_LIBRARIES};${METIS_LIBRARIES};${LIB_DOWNLOAD_PATH}/lib/superlu_dist.lib;${LAPACK_LIBRARIES};${BLAS_LIBRARIES};${PARMETIS_LIBRARIES}" CACHE PATH "SUPERLU_LIBRARIES is set" FORCE)
					else( WIN32 )
						include(GNUInstallDirs)
						message("gnu install library dirs ${CMAKE_INSTALL_LIBDIR}")
						set(SUPERLU_LIBRARIES "${SUPERLU_LIBRARIES};${METIS_LIBRARIES};${LIB_DOWNLOAD_PATH}/${CMAKE_INSTALL_LIBDIR}/libsuperlu_dist.a;${LAPACK_LIBRARIES};${BLAS_LIBRARIES};${PARMETIS_LIBRARIES}" CACHE PATH "SUPERLU_LIBRARIES is set" FORCE)
					endif( WIN32 )
					message("-- SUPERLU_INCLUDE_DIR: ${SUPERLU_INCLUDES}")
					message("-- SUPERLU_LIBRARIES:   ${SUPERLU_LIBRARIES}")					
				else()
					message("-- Skip SuperLU_DIST as either USE_MPI or USE_PARTITIONER_PARMETIS are not checked for compilation")
				endif()
				
				
			else(USE_SOLVER_SUPERLU_DOWNLOAD)
				set(USE_SOLVER_SUPERLU OFF CACHE BOOL "Use SuperLU solver" FORCE)
			endif(USE_SOLVER_SUPERLU_DOWNLOAD)
		else()
			MESSAGE("\nFound SuperLU_DIST!  Here are the details: ")
			set(SUPERLU_INCLUDES "${SUPERLU_INCLUDES};${SUPERLU_DIST_INCLUDES}")
			set(SUPERLU_LIBRARIES "${SUPERLU_LIBRARIES};${SUPERLU_DIST_LIBRARIES}")
			MESSAGE("INCLUDES:  ${SUPERLU_INCLUDES}")
			MESSAGE("LIBRARIES: ${SUPERLU_LIBRARIES}")
		endif()
	endif()
	if( USE_SOLVER_SUPERLU )
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
		Source/Headers/inmost_optimizer.h
                           Source/Headers/container.hpp)

#if( COMPILE_EXAMPLES )
#list(APPEND INMOST_INSTALL_HEADERS "${PROJECT_SOURCE_DIR}/Examples/AdaptiveMesh/amesh.h")
#message(${INMOST_INSTALL_HEADERS})
#endif( COMPILE_EXAMPLES )

include(CPack)
export(TARGETS inmost FILE inmost-targets.cmake)
export(PACKAGE inmost)
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/Source/Headers" "${PROJECT_BINARY_DIR}")
set(CONF_LIBRARY_DIRS "")
if( COMPILE_EXAMPLES )
list(APPEND CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/Examples/AdaptiveMesh")
list(APPEND CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/Examples/GridTools")
list(APPEND CONF_LIBRARY_DIRS "${PROJECT_BINARY_DIR}/Examples/AdaptiveMesh")
list(APPEND CONF_LIBRARY_DIRS "${PROJECT_BINARY_DIR}/Examples/GridTools")
endif( COMPILE_EXAMPLES )
configure_file(inmost-config.cmake.in "${PROJECT_BINARY_DIR}/inmost-config.cmake" @ONLY)
set(CONF_INCLUDE_DIRS "\${inmost_DIR}/include")
set(CONF_LIBRARY_DIRS "\$")
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
	"${PROJECT_SOURCE_DIR}/Source/Headers/inmost_optimizer.h"
	"${PROJECT_SOURCE_DIR}/Source/Headers/container.hpp")

#if( COMPILE_EXAMPLES )
#set_property(TARGET inmost APPEND PROPERTY PUBLIC_HEADER "${PROJECT_SOURCE_DIR}/Examples/AdaptiveMesh/amesh.h")
#get_property(PRINT_PUBLIC_HEADER TARGET inmost PROPERTY PUBLIC_HEADER)
#message(${PRINT_PUBLIC_HEADER})
#endif( COMPILE_EXAMPLES )

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
