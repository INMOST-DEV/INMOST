project(OctreeCutcell)
set(SOURCE main.cpp
           rotate.cpp
           rotate.h
	   tree.cpp
	   tree.h
	   triang.cpp
	   triang.h
	   gettype.hpp
	   obj.cpp
	   obj.h
	   octgrid.cpp
	   octgrid.h
	   point.h
	   proj.hpp
           my_glut.h)


add_executable(OctreeCutcell ${SOURCE})

target_link_libraries(OctreeCutcell inmost)


if(USE_MPI)
  message("linking with mpi")
  target_link_libraries(OctreeCutcell ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    message("adding link flags")
    set_target_properties(OctreeCutcell PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)

find_package(OpenGL)
find_package(GLUT)

if(OPENGL_FOUND)
	if(GLUT_FOUND)
		include_directories(${OPENGL_INCLUDE_DIR})
		include_directories(${GLUT_INCLUDE_DIR})
		add_definitions(-D__GRAPHICS__)
		target_link_libraries(OctreeCutcell ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
	else(GLUT_FOUND)
		message("GLUT not found")
	endif(GLUT_FOUND)
else(OPENGL_FOUND)
	message("OpenGL not found")
endif(OPENGL_FOUND)

