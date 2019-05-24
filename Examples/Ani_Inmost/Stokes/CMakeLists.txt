#cmake_minimum_required(VERSION 2.6)
#enable_language(Fortran)
project(Stokes)
#enable_language (CXX C Fortran)


#set(CMAKE_CXX_STANDARD 11)


#set(CMAKE_FORTRAN_FLAGS "${CMAKE_FORTRAN_FLAGS} -O3")
#set(CMAKE_Fortran_FLAGS_DEBUG  "${CMAKE_FORTRAN_FLAGS_DEBUG} -O0 -g")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
#set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0  -g")
#set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0  -g")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")



#set(INMOST_SRC "/home/vasily/Work/INMOST_ANI/INMOST_MAIN")
#set(INMOST_BUILD "/home/vasily/Work/INMOST_ANI/build_main")
#set(PETSC_DIR "/home/vasily/Work/INMOST_ANI/new_generation/petsc-3.5.4")
#set(PETSC_ARCH "arch-linux2-c-opt")
#set (ANI_DIR "/home/vasily/Work/INMOST_ANI/new_generation/ani3D-3.1")
#set( PARENT_DIR "../..")
#include_directories(${INMOST_SRC}/Source/Headers)
#include_directories(${INMOST_BUILD})


#include_directories(${PETSC_DIR}/include)
#include_directories(${PETSC_DIR}/${PETSC_ARCH}/include)
include_directories($ENV{ANI_DIR}/src/aniFEM)

#include_directories(../../)

#add_library (fort_lib STATIC  forlibfem.f)
#SET_TARGET_PROPERTIES(fort_lib PROPERTIES LINKER_LANGUAGE Fortran)





set(EXEC_NAME Stokes)

add_executable(${EXEC_NAME} main.cpp forlibfem.f)





#target_link_libraries(${EXEC_NAME}  fort_lib)
target_link_libraries(${EXEC_NAME} inmost_ani)
target_link_libraries(${EXEC_NAME} refine)
target_link_libraries(${EXEC_NAME} inmost)


#link_directories(${ANI_DIR}/lib)
#link_directories(${PETSC_DIR}/${PETSC_ARCH}/lib/)

#target_link_libraries(${EXEC_NAME} mba3D-3.0)
#target_link_libraries(${EXEC_NAME} fem3D-3.0)
target_link_libraries(${EXEC_NAME} $ENV{ANI_DIR}/lib/libmba3D-3.0.a)
target_link_libraries(${EXEC_NAME} $ENV{ANI_DIR}/lib/libfem3D-3.0.a)
if(USE_SOLVER_PETSC)
else()
target_link_libraries(${EXEC_NAME} $ENV{ANI_DIR}/lib/liblapack-3.2.a)
target_link_libraries(${EXEC_NAME} $ENV{ANI_DIR}/lib/liblapack_ext-3.2.a)
target_link_libraries(${EXEC_NAME} $ENV{ANI_DIR}/lib/libblas-3.2.a)
endif()

if(USE_MPI)
  message("linking ${EXEC_NAME} with MPI")
  target_link_libraries(${EXEC_NAME} ${MPI_LIBRARIES})
  if(MPI_LINK_FLAGS)
    set_target_properties(${EXEC_NAME} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()
endif(USE_MPI)

if(USE_SOLVER)
  if(USE_SOLVER_ANI)
    message("linking ${EXEC_NAME} with ani3d and BLAS")
    target_link_libraries(${EXEC_NAME} ani3d ${BLAS_LIBRARIES})
    if(BLAS_LINKER_FLAGS)
      set_target_properties(${EXEC_NAME} PROPERTIES LINK_FLAGS "${BLAS_LINKER_FLAGS}")
    endif()
  endif()
  if(USE_SOLVER_PETSC)
    message("linking ${EXEC_NAME} with PETSc")
    target_link_libraries(${EXEC_NAME} ${PETSC_LIBRARIES})
  endif()
  if(USE_SOLVER_TRILINOS)
    message("linking ${EXEC_NAME} with Trilinos")
    target_link_libraries(${EXEC_NAME} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
  endif()
  if(USE_SOLVER_METIS)
    message("linking ${EXEC_NAME} with Metis")
    target_link_libraries(${EXEC_NAME} ${METIS_LIBRARIES})
  endif()
  if(USE_SOLVER_MONDRIAAN)
    message("linking ${EXEC_NAME} with Mondriaan")
    target_link_libraries(${EXEC_NAME} ${MONDRIAAN_LIBRARIES})
  endif()
  if(USE_SOLVER_SUPERLU)
    message("linking ${EXEC_NAME} with SuperLU")
    target_link_libraries(${EXEC_NAME} ${SUPERLU_LIBRARIES})
  endif()
endif()


if(USE_PARTITIONER)
  if(USE_PARTITIONER_ZOLTAN)
    message("linking ${EXEC_NAME} with Zoltan")
    target_link_libraries(${EXEC_NAME} ${ZOLTAN_LIBRARIES})
  endif()
  if(USE_PARTITIONER_PARMETIS)
    message("linking ${EXEC_NAME} with ParMETIS")
    target_link_libraries(${EXEC_NAME} ${PARMETIS_LIBRARIES})
  endif()
endif()


#target_link_libraries(${EXEC_NAME} ${INMOST_BUILD}/libinmost.a)
#message("${PETSC_DIR}/${PETSC_ARCH}")
#target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.so)

#target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libparmetis.so)
#target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libmetis.so)
#target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libf2clapack.a)
#target_link_libraries(${EXEC_NAME} ${PETSC_DIR}/${PETSC_ARCH}/lib/libf2cblas.a)


#target_link_libraries(${EXEC_NAME} petsc)

#target_link_libraries(${EXEC_NAME} parmetis)
#target_link_libraries(${EXEC_NAME} metis)
#target_link_libraries(${EXEC_NAME} f2clapack)
#target_link_libraries(${EXEC_NAME} f2cblas)

#target_link_libraries(${EXEC_NAME} -ldl)
target_link_libraries(${EXEC_NAME} -lgfortran)
#target_link_libraries(${EXEC_NAME} -lX11)
#target_link_libraries(${EXEC_NAME} -lssl)
#target_link_libraries(${EXEC_NAME} -lcrypto)
