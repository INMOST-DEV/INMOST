add_executable(mpi_window_test000 main.cpp)
target_link_libraries(mpi_window_test000 inmost ${MPI_LIBRARIES})
if(MPI_LINK_FLAGS)
  set_target_properties(mpi_window_test000 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

if(MPIEXEC_EXECUTABLE)
  foreach(ranks 2 4)
    add_test(NAME mpi_window_test000_np_${ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
              ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mpi_window_test000> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(mpi_window_test000_np_${ranks} PROPERTIES
      TIMEOUT 60 PROCESSORS ${ranks})
  endforeach()
endif()
