project(autodiff_test000)
set(SOURCE main.cpp)

add_executable(container_test000 ${SOURCE})
target_link_libraries(container_test000 inmost)

if(USE_MPI)
  message("linking container_test000 with MPI")
  target_link_libraries(container_test000 ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(container_test000 PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)

if(USE_ZLIB)
  message("linking container_test000 with ZLIB")
  target_link_libraries(container_test000 ${ZLIB_LIBRARIES}) 
endif(USE_ZLIB)


add_test(NAME container_test000_pool_array COMMAND $<TARGET_FILE:container_test000> 0)

