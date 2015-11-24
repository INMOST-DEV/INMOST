project(pmesh_test000)
set(SOURCE main.cpp)

add_executable(pmesh_test000 ${SOURCE})
target_link_libraries(pmesh_test000 inmost)

if(USE_MPI)
  message("linking pmesh_test000 with MPI")
  target_link_libraries(pmesh_test000 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(pmesh_test000 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)

if( USE_MPI AND EXISTS ${MPIEXEC} )
  add_test(NAME pmesh_test000_parallel_np_2  COMMAND ${MPIEXEC} -np 2 $<TARGET_FILE:pmesh_test000> ${CMAKE_CURRENT_SOURCE_DIR}/3x3x3.vtk)
  add_test(NAME pmesh_test000_parallel_np_3  COMMAND ${MPIEXEC} -np 3 $<TARGET_FILE:pmesh_test000> ${CMAKE_CURRENT_SOURCE_DIR}/3x3x3.vtk)
  add_test(NAME pmesh_test000_parallel_np_4  COMMAND ${MPIEXEC} -np 4 $<TARGET_FILE:pmesh_test000> ${CMAKE_CURRENT_SOURCE_DIR}/3x3x3.vtk)
endif()
