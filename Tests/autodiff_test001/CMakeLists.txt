project(autodiff_test001)
set(SOURCE main.cpp)

add_executable(autodiff_test001 ${SOURCE})
target_link_libraries(autodiff_test001 inmost)

if(USE_MPI)
  message("linking autodiff_test001 with MPI")
  target_link_libraries(autodiff_test001 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(autodiff_test001 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)



add_test(NAME autodiff_test001_hessian_0 COMMAND $<TARGET_FILE:autodiff_test001> 0)
