project(AdaptiveMesh)

set(LIBSOURCE amesh.cpp amesh.h)
set(SOURCE main.cpp)

add_library(AdaptiveMeshLib ${LIBSOURCE})
add_executable(AdaptiveMesh ${SOURCE})

target_link_libraries(AdaptiveMesh inmost AdaptiveMeshLib)

if(USE_MPI)
  message("linking AdaptiveMesh with MPI")
  target_link_libraries(AdaptiveMesh ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    set_target_properties(AdaptiveMesh PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)

set_property(TARGET AdaptiveMeshLib PROPERTY PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/amesh.h")

install(TARGETS AdaptiveMeshLib EXPORT inmost-targets 
        LIBRARY DESTINATION lib 
        ARCHIVE DESTINATION lib
        PUBLIC_HEADER DESTINATION include)
install(TARGETS AdaptiveMesh EXPORT inmost-targets RUNTIME DESTINATION bin)