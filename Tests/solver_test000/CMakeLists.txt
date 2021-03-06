project(solver_test000)
set(SOURCE main.cpp)

add_executable(solver_test000 ${SOURCE})
target_link_libraries(solver_test000 inmost)

if(USE_MPI)
  message("linking solver_test000 with MPI")
  target_link_libraries(solver_test000 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(solver_test000 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)


if(USE_SOLVER_ANI)
  message("linking solver_test000 with ani3d and BLAS")
  target_link_libraries(solver_test000 ani3d ${BLAS_LIBRARIES})
  if(BLAS_LINKER_FLAGS)
    set_target_properties(solver_test000 PROPERTIES LINK_FLAGS "${BLAS_LINKER_FLAGS}")
  endif()
endif()
if(USE_SOLVER_PETSC)
  message("linking solver_test000 with PETSc")
  target_link_libraries(solver_test000 ${PETSC_LIBRARIES})
endif()
if(USE_SOLVER_TRILINOS AND USE_MPI)
  message("linking solver_test000 with Trilinos")
  target_link_libraries(solver_test000 ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
endif()
if(USE_SOLVER_METIS)
  message("linking solver_test000 with Metis")
  target_link_libraries(solver_test000 ${METIS_LIBRARIES})
endif()
if(USE_SOLVER_MONDRIAAN)
  message("linking solver_test000 with Mondriaan")
  target_link_libraries(solver_test000 ${MONDRIAAN_LIBRARIES})
endif()
if(USE_SOLVER_SUPERLU)
  message("linking solver_test000 with SuperLU")
  target_link_libraries(solver_test000 ${SUPERLU_LIBRARIES})
endif()


add_test(NAME solver_test000_serial_inner_ilu2                 COMMAND $<TARGET_FILE:solver_test000> 0 inner_ilu2)
add_test(NAME solver_test000_serial_inner_ddpqiluc             COMMAND $<TARGET_FILE:solver_test000> 0 inner_ddpqiluc2)

add_test(NAME solver_test000_serial_inner_mptiluc             COMMAND $<TARGET_FILE:solver_test000> 0 inner_mptiluc)
add_test(NAME solver_test000_serial_inner_mlmptiluc             COMMAND $<TARGET_FILE:solver_test000> 0 inner_mlmptiluc)

add_test(NAME solver_test000_serial_inner_mptilu2             COMMAND $<TARGET_FILE:solver_test000> 0 inner_mptilu2)


if(USE_SOLVER_PETSC)
add_test(NAME solver_test000_serial_petsc                      COMMAND $<TARGET_FILE:solver_test000> 0 petsc)
add_test(NAME solver_test000_copy_serial_petsc                      COMMAND $<TARGET_FILE:solver_test000> 0 petsc 1)
add_test(NAME solver_test000_assign_serial_petsc                      COMMAND $<TARGET_FILE:solver_test000> 0 petsc 2)
endif()

if(USE_SOLVER_TRILINOS AND USE_MPI)
add_test(NAME solver_test000_serial_trilinos_aztec             COMMAND $<TARGET_FILE:solver_test000> 0 trilinos_aztec)
add_test(NAME solver_test000_serial_trilinos_ifpack            COMMAND $<TARGET_FILE:solver_test000> 0 trilinos_ifpack)
add_test(NAME solver_test000_serial_trilinos_ml                COMMAND $<TARGET_FILE:solver_test000> 0 trilinos_ml)
add_test(NAME solver_test000_serial_trilinos_belos             COMMAND $<TARGET_FILE:solver_test000> 0 trilinos_belos)
endif()

if(USE_SOLVER_ANI)
add_test(NAME solver_test000_serial_ani COMMAND $<TARGET_FILE:solver_test000> 0 ani)
endif()

if(SOLVER_DEFINITIONS MATCHES "^.*HAVE_SOLVER_FCBIILU2.*$")
  add_test(NAME solver_test000_serial_fcbiilu2 COMMAND $<TARGET_FILE:solver_test000> 0 fcbiilu2)
endif()

if( USE_MPI )
if( EXISTS ${MPIEXEC} )
add_test(NAME solver_test000_parallel_normal_inner_ilu2        COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 0 inner_ilu2)
add_test(NAME solver_test000_parallel_permute1_inner_ilu2      COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 1 inner_ilu2)
add_test(NAME solver_test000_parallel_permute2_inner_ilu2      COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 2 inner_ilu2)

add_test(NAME solver_test000_parallel_normal_inner_ddpqiluc      COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 0 inner_ddpqiluc2)
add_test(NAME solver_test000_parallel_permute1_inner_ddpqiluc    COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 1 inner_ddpqiluc2)
add_test(NAME solver_test000_parallel_permute2_inner_ddpqiluc    COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 2 inner_ddpqiluc2)

add_test(NAME solver_test000_parallel_normal_inner_mptiluc      COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 0 inner_mptiluc)
add_test(NAME solver_test000_parallel_permute1_inner_mptiluc    COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 1 inner_mptiluc)
add_test(NAME solver_test000_parallel_permute2_inner_mptiluc    COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 2 inner_mptiluc)

add_test(NAME solver_test000_parallel_normal_inner_mlmptiluc      COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 0 inner_mlmptiluc)
add_test(NAME solver_test000_parallel_permute1_inner_mlmptiluc    COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 1 inner_mlmptiluc)
add_test(NAME solver_test000_parallel_permute2_inner_mlmptiluc    COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 2 inner_mlmptiluc)

add_test(NAME solver_test000_parallel_normal_inner_mptilu2      COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 0 inner_mptilu2)
add_test(NAME solver_test000_parallel_permute1_inner_mptilu2    COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 1 inner_mptilu2)
add_test(NAME solver_test000_parallel_permute2_inner_mptilu2    COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 2 inner_mptilu2)

if(USE_SOLVER_PETSC)
add_test(NAME solver_test000_parallel_normal_petsc             COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 0 petsc)
add_test(NAME solver_test000_copy_parallel_normal_petsc             COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 0 petsc 1)
add_test(NAME solver_test000_assign_parallel_normal_petsc             COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 0 petsc 2)
#add_test(NAME solver_test000_parallel_permute1_petsc           COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 1 2)
#add_test(NAME solver_test000_parallel_permute2_petsc           COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 2 2)
endif()

if(USE_SOLVER_TRILINOS)
add_test(NAME solver_test000_parallel_normal_trilinos_aztec    COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 0 trilinos_aztec)
add_test(NAME solver_test000_parallel_permute1_trilinos_aztec  COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 1 trilinos_aztec)
add_test(NAME solver_test000_parallel_permute2_trilinos_aztec  COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 2 trilinos_aztec)

add_test(NAME solver_test000_parallel_normal_trilinos_ifpack   COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 0 trilinos_ifpack)
add_test(NAME solver_test000_parallel_permute1_trilinos_ifpack COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 1 trilinos_ifpack)
add_test(NAME solver_test000_parallel_permute2_trilinos_ifpack COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 2 trilinos_ifpack)

add_test(NAME solver_test000_parallel_normal_trilinos_ml       COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 0 trilinos_ml)
add_test(NAME solver_test000_parallel_permute1_trilinos_ml     COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 1 trilinos_ml)
add_test(NAME solver_test000_parallel_permute2_trilinos_ml     COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 2 trilinos_ml)

#add_test(NAME solver_test000_parallel_normal_trilinos_belos    COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 0 6)
#add_test(NAME solver_test000_parallel_permute1_trilinos_belos  COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 1 6)
#add_test(NAME solver_test000_parallel_permute2_trilinos_belos  COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 2 6)
endif()

if(SOLVER_DEFINITIONS MATCHES "^.*HAVE_SOLVER_FCBIILU2.*$")
  add_test(NAME solver_test000_parallel_normal_fcbiilu2             COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 0 fcbiilu2)
endif()

if(SOLVER_DEFINITIONS MATCHES "^.*HAVE_SOLVER_K3BIILU2.*$")
  add_test(NAME solver_test000_parallel_normal_k3biilu2             COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:solver_test000> 0 k3biilu2)
endif()

endif()
endif()
