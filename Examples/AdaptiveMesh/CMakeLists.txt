project(AdaptiveMesh)


add_library(AdaptiveMeshLib amesh.cpp amesh.h)
add_executable(AdaptiveMesh main.cpp)
add_executable(AdaptiveMeshCube main_cube.cpp)
add_executable(AdaptiveMeshSphere main_sphere.cpp)
#add_executable(AdaptiveMeshTestSet test_set.cpp)

target_link_libraries(AdaptiveMesh inmost AdaptiveMeshLib)
target_link_libraries(AdaptiveMeshCube inmost AdaptiveMeshLib)
target_link_libraries(AdaptiveMeshSphere inmost AdaptiveMeshLib)
#target_link_libraries(AdaptiveMeshTestSet inmost AdaptiveMeshLib)

if(USE_MPI)
  message("linking AdaptiveMesh with MPI")
  target_link_libraries(AdaptiveMesh ${MPI_LIBRARIES}) 
  target_link_libraries(AdaptiveMeshCube ${MPI_LIBRARIES}) 
  target_link_libraries(AdaptiveMeshSphere ${MPI_LIBRARIES}) 
  #target_link_libraries(AdaptiveMeshTestSet ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(AdaptiveMesh PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    set_target_properties(AdaptiveMeshCube PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    set_target_properties(AdaptiveMeshSphere PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    #set_target_properties(AdaptiveMeshTestSet PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)


if(USE_PARTITIONER)
  if(USE_PARTITIONER_ZOLTAN)
    message("linking AdaptiveMesh with Zoltan")
    target_link_libraries(AdaptiveMesh ${ZOLTAN_LIBRARIES})
    target_link_libraries(AdaptiveMeshCube ${ZOLTAN_LIBRARIES}) 
    target_link_libraries(AdaptiveMeshSphere ${ZOLTAN_LIBRARIES}) 
  endif()
  if(USE_PARTITIONER_PARMETIS)
    message("linking AdaptiveMesh with ParMETIS")
    target_link_libraries(AdaptiveMesh ${PARMETIS_LIBRARIES})
    target_link_libraries(AdaptiveMeshCube ${PARMETIS_LIBRARIES})
    target_link_libraries(AdaptiveMeshSphere ${PARMETIS_LIBRARIES})
  endif()
endif()

set_property(TARGET AdaptiveMeshLib PROPERTY PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/amesh.h")

install(TARGETS AdaptiveMeshLib EXPORT inmost-targets 
        LIBRARY DESTINATION lib 
        ARCHIVE DESTINATION lib
        PUBLIC_HEADER DESTINATION include)
install(TARGETS AdaptiveMesh EXPORT inmost-targets RUNTIME DESTINATION bin)
install(TARGETS AdaptiveMeshCube EXPORT inmost-targets RUNTIME DESTINATION bin)
install(TARGETS AdaptiveMeshSphere EXPORT inmost-targets RUNTIME DESTINATION bin)
#install(TARGETS AdaptiveMeshTestSet EXPORT inmost-targets RUNTIME DESTINATION bin)
