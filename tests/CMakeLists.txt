if(USE_SOLVER)
add_subdirectory(solver_test000)
endif(USE_SOLVER)

if(USE_MESH AND USE_MPI)
add_subdirectory(pmesh_test000)
endif()