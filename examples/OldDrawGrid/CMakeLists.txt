project(OldDrawGrid)
set(SOURCE main.cpp
           rotate.cpp
           rotate.h
           my_glut.h)
find_package(OpenGL)

if (WIN32)
  find_path(GLUT_INCLUDE_DIR GL/glut.h )
  if (CYGWIN)
    find_library(GLUT_glut_LIBRARY glut32 )
  else ()
    if(BORLAND)
      set (GLUT_glut_LIBRARY import32 CACHE STRING "GLUT library for win32")
    else()
      set (GLUT_glut_LIBRARY opengl32 CACHE STRING "GLUT library for win32")
    endif()
  endif ()
  if(NOT GLUT_glut_LIBRARY)
    #-- use old way to find GLUT
    find_path( GLUT_INCLUDE_DIR NAMES GL/glut.h
      PATHS  ${GLUT_ROOT_PATH}/include )
    find_library( GLUT_glut_LIBRARY NAMES glut glut32 freeglut
      PATHS
      ${OPENGL_LIBRARY_DIR}
      ${GLUT_ROOT_PATH}/Release
      )
  endif(NOT GLUT_glut_LIBRARY)

  #include(FindPackageHandleStandardArgs.cmake)
  #FIND_PACKAGE_HANDLE_STANDARD_ARGS(GLUT REQUIRED_VARS GLUT_glut_LIBRARY GLUT_INCLUDE_DIR)

  if (GLUT_glut_LIBRARY)
    message("GLUT library found")
    set( GLUT_LIBRARIES ${GLUT_glut_LIBRARY})
    #The following deprecated settings are for backwards compatibility with CMake1.4
    set (GLUT_LIBRARY ${GLUT_LIBRARIES})
	if(NOT GLUT_INCLUDE_DIR)
		message("GLUT include directory not found - trying to set OpenGL library")
		message("OpenGL library is ${OPENGL_INCLUDE_DIR}")
		set(GLUT_INCLUDE_DIR ${OPENGL_INCLUDE_DIR})
        	set(GLUT_INCLUDE_PATH ${GLUT_INCLUDE_DIR})
		set(GLUT_FOUND 1)
	else(NOT GLUT_INCLUDE_DIR)
		set(GLUT_FOUND 1)
	endif(NOT GLUT_INCLUDE_DIR)
  endif(GLUT_glut_LIBRARY)

else(WIN32)
  find_package(GLUT)
endif(WIN32)

if(OPENGL_FOUND)
	include_directories(${OpenGL_INCLUDE_DIRS})
	link_directories(${OpenGL_LIBRARY_DIRS})
	add_definitions(${OpenGL_DEFINITIONS})
	if(GLUT_FOUND)
		include_directories(${GLUT_INCLUDE_DIRS})
		link_directories(${GLUT_LIBRARY_DIRS})
		add_definitions(${GLUT_DEFINITIONS})
		add_executable(OldDrawGrid ${SOURCE})
		target_link_libraries(OldDrawGrid ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} inmost)
		if(USE_MPI)
		  target_link_libraries(OldDrawGrid ${MPI_LIBRARIES}) 
		  if(MPI_LINK_FLAGS)
            set_target_properties(OldDrawGrid PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
            endif() 
		endif(USE_MPI)
	else(GLUT_FOUND)
		message("GLUT not found")
	endif(GLUT_FOUND)
else(OPENGL_FOUND)
	message("OpenGL not found")
endif(OPENGL_FOUND)

