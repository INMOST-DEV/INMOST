project(DrawGrid)
set(SOURCE main.cpp
           rotate.cpp
           clipboard.cpp
           rotate.h
           inc_glut.h
           clipboard.h
           color.h
           coord.h
           coord.cpp
           octree.h
           streamline.h
           streamline.cpp
           svg_line.h
           svg_line.cpp
           face2gl.h
           face2gl.cpp
		   color_bar.h
		   color_bar.cpp
		   printtext.cpp
		   tga.h
		   tga.cpp
		   screenshot.h
		   screenshot.cpp
		   volumetric.h
		   volumetric.cpp
		   input.h
		   input.cpp
		   picker.h
		   picker.cpp
		   clipper.h
		   clipper.cpp
		   vector.h
		   vector.cpp)
find_package(OpenGL)
find_package(GLUT)

if(OPENGL_FOUND)
  if(GLUT_FOUND)
    message("linking DrawGrid with GLUT and OpenGL")
    include_directories(${OPENGL_INCLUDE_DIR})
    include_directories(${GLUT_INCLUDE_DIR})
    add_executable(DrawGrid ${SOURCE})
    target_link_libraries(DrawGrid inmost ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
    if(USE_MPI)
      target_link_libraries(DrawGrid ${MPI_LIBRARIES}) 
      if(MPI_LINK_FLAGS)
        set_target_properties(DrawGrid PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
      endif() 
    endif(USE_MPI)
    install(TARGETS DrawGrid EXPORT inmost-targets RUNTIME DESTINATION bin)
  else(GLUT_FOUND)
    message("GLUT not found, not building DrawGrid")
  endif(GLUT_FOUND)
else(OPENGL_FOUND)
  message("OpenGL not found, not building DrawGrid")
endif(OPENGL_FOUND)
