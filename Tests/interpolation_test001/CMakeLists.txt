project(interpolation_test001)
set(SOURCE main.cpp)

add_executable(interpolation_test001 ${SOURCE})
target_link_libraries(interpolation_test001 inmost)

if(USE_MPI)
  message("linking interpolation_test001 with MPI")
  target_link_libraries(interpolation_test001 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(interpolation_test001 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)

if(USE_ZLIB)
  message("linking interpolation_test001 with ZLIB")
  target_link_libraries(interpolation_test001 ${ZLIB_LIBRARIES}) 
endif(USE_ZLIB)

add_test(NAME interpolation_test001_cell         COMMAND $<TARGET_FILE:interpolation_test001> 0)
