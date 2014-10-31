project(GridGen)
set(SOURCE main.cpp)

add_executable(GridGen ${SOURCE})

target_link_libraries(GridGen inmost)

if(USE_MPI)
  message("linking GridGen with MPI")
  target_link_libraries(GridGen ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(GridGen PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
