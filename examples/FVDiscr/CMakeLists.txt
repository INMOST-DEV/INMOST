project(FVDiscr)
set(SOURCE main.cpp)

add_executable(FVDiscr ${SOURCE})

target_link_libraries(FVDiscr inmost)

if(USE_SOLVER_ANI)
  message("linking FVDiscr with ani3d")
  target_link_libraries(FVDiscr ani3d)
endif()

if(USE_SOLVER_PETSC)
  message("linking FVDiscr with PETSc")
  link_directories(${PETSC_LIBRARY_DIRS})
  target_link_libraries(FVDiscr ${PETSC_LIBRARIES})
endif()

if(USE_PARTITIONER)
  if(USE_PARTITIONER_ZOLTAN)
    target_link_libraries(FVDiscr ${ZOLTAN_LIBRARIES})
  endif()
  if(USE_PARTITIONER_PARMETIS)
    target_link_libraries(FVDiscr ${PARMETIS_LIBRARIES})
  endif()
endif()

if(USE_MPI)
  message("linking FVDiscr with MPI")
  target_link_libraries(FVDiscr ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(FVDiscr PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)
