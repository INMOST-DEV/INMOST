project(DrawGrid)
set(SOURCE main.cpp
           rotate.cpp
           rotate.h
           my_glut.h)
find_package(OpenGL)
find_package(GLUT)

if(OPENGL_FOUND)
	if(GLUT_FOUND)
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
	else(GLUT_FOUND)
		message("GLUT not found")
	endif(GLUT_FOUND)
else(OPENGL_FOUND)
	message("OpenGL not found")
endif(OPENGL_FOUND)

install(TARGETS DrawGrid EXPORT inmost-targets RUNTIME DESTINATION bin)