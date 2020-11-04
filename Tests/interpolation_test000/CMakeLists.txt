project(interpolation_test000)
set(SOURCE main.cpp)

add_executable(interpolation_test000 ${SOURCE})
target_link_libraries(interpolation_test000 inmost)

if(USE_MPI)
  message("linking interpolation_test000 with MPI")
  target_link_libraries(interpolation_test000 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(interpolation_test000 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)


add_test(NAME interpolation_test000_face         COMMAND $<TARGET_FILE:interpolation_test000> 0)
