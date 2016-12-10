if(USE_AUTODIFF)
add_subdirectory(autodiff_test000)
endif(USE_AUTODIFF)

if(USE_SOLVER)
add_subdirectory(solver_test000)
#add_subdirectory(solver_test001)
add_subdirectory(solver_test002)
endif(USE_SOLVER)

if(USE_MESH AND USE_MPI)
add_subdirectory(pmesh_test000)
if(USE_PARTITIONER)
add_subdirectory(pmesh_test001)
endif()
endif()

add_subdirectory(xml_reader_test000)
