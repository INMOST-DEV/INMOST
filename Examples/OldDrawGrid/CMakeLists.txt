project(OldDrawGrid)
set(SOURCE main.cpp
           rotate.cpp
           clipboard.cpp
           rotate.h
           my_glut.h
           clipboard.h)
find_package(OpenGL)
find_package(GLUT)

if(OPENGL_FOUND)
  if(GLUT_FOUND)
    message("linking OldDrawGrid with GLUT and OpenGL")
    include_directories(${OPENGL_INCLUDE_DIR})
    include_directories(${GLUT_INCLUDE_DIR})
    add_executable(OldDrawGrid ${SOURCE})
    target_link_libraries(OldDrawGrid inmost ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
    if(USE_MPI)
      target_link_libraries(OldDrawGrid ${MPI_LIBRARIES})
      if(MPI_LINK_FLAGS)
        set_target_properties(OldDrawGrid PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
      endif()
    endif(USE_MPI)
    install(TARGETS OldDrawGrid EXPORT inmost-targets RUNTIME DESTINATION bin)
  else(GLUT_FOUND)
    message("GLUT not found, not building OldDrawGrid")
  endif(GLUT_FOUND)
else(OPENGL_FOUND)
  message("OpenGL not found, not building OldDrawGrid")
endif(OPENGL_FOUND)
