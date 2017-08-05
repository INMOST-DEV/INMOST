project(geom_test000)
set(SOURCE main.cpp)

add_executable(geom_test000 ${SOURCE})
target_link_libraries(geom_test000 inmost)

if(USE_MPI)
  message("linking geom_test000 with MPI")
  target_link_libraries(geom_test000 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(geom_test000 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)

add_test(NAME geom_test000_dual4  COMMAND $<TARGET_FILE:geom_test000> ${CMAKE_CURRENT_SOURCE_DIR}/d4.pmf)
add_test(NAME geom_test000_cube4  COMMAND $<TARGET_FILE:geom_test000> ${CMAKE_CURRENT_SOURCE_DIR}/c4.pmf)
add_test(NAME geom_test000_tetra4  COMMAND $<TARGET_FILE:geom_test000> ${CMAKE_CURRENT_SOURCE_DIR}/t4.pmf)
add_test(NAME geom_test000_dual4_split  COMMAND $<TARGET_FILE:geom_test000> ${CMAKE_CURRENT_SOURCE_DIR}/d4sp.pmf)
