#cmake_minimum_required(VERSION 2.6)
enable_language(Fortran)
#project(INMOST_ANI_LIB)
enable_language (CXX C Fortran)

#set(CMAKE_CXX_STANDARD 11)


#set(CMAKE_FORTRAN_FLAGS "${CMAKE_FORTRAN_FLAGS} -O3")
#set(CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_FORTRAN_FLAGS_DEBUG} -O0 -g")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0  -g")
#set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0  -g")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#message(${CMAKE_C_FLAGS})

#set(INMOST_SRC "/home/vasily/Work/INMOST_ANI/new_generation/INMOST")
#set(INMOST_BUILD "/home/vasily/Work/INMOST_ANI/new_generation/build")
#option(USE_MPI "Compile with MPI support" ON)

#include_directories(${INMOST_SRC}/Source/Headers)
#include_directories(${INMOST_BUILD})

add_library (refine STATIC ./inmost_ani_lib/refine.f)
SET_TARGET_PROPERTIES(refine PROPERTIES LINKER_LANGUAGE Fortran)

list(APPEND SOURCE_FILES ./inmost_ani_lib/inmost_ani_mesh.cpp)
list(APPEND SOURCE_FILES ./inmost_ani_lib/inmost_ani_mesh.h)
list(APPEND SOURCE_FILES ./inmost_ani_lib/inmost_ani_fem.cpp)
list(APPEND SOURCE_FILES ./inmost_ani_lib/inmost_ani_fem.h)
list(APPEND SOURCE_FILES ./inmost_ani_lib/utils.cpp)
list(APPEND SOURCE_FILES ./inmost_ani_lib/utils.h)
add_library (inmost_ani STATIC ${SOURCE_FILES})
include_directories(./inmost_ani_lib)

if(USE_MPI)
#	find_package(MPI)
	if(NOT MPI_FOUND)
#		set(USE_MPI OFF CACHE BOOL "Compile with MPI support" FORCE)
		message("MPI NOT FOUND")
  	else()
		include_directories(${MPI_INCLUDE_PATH})
#		set_target_properties(inmost_ani PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
#		message("MPI FOUND")
	endif()
endif()
add_subdirectory(Unst_CD_Proj)
add_subdirectory(Diffusion_P2)
add_subdirectory(Elasticity)
add_subdirectory(Stokes)
add_subdirectory(MixedHybrid)


