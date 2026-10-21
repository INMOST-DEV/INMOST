add_executable(pmf_test000 main.cpp)
target_link_libraries(pmf_test000 inmost ${MPI_LIBRARIES})
set_property(TARGET pmf_test000 PROPERTY CXX_STANDARD 11)
if(MPI_LINK_FLAGS)
  set_target_properties(pmf_test000 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
if(MPIEXEC_EXECUTABLE)
  foreach(ranks 2 4)
    add_test(NAME pmf_test000_np_${ranks}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
              ${MPIEXEC_PREFLAGS} $<TARGET_FILE:pmf_test000>
              ${CMAKE_CURRENT_BINARY_DIR}/roundtrip-${ranks} ${MPIEXEC_POSTFLAGS})
    set_tests_properties(pmf_test000_np_${ranks} PROPERTIES TIMEOUT 120 PROCESSORS ${ranks})
  endforeach()
endif()
