project(Octree)
set(SOURCE main.cpp
           rotate.cpp
           rotate.h
	   octgrid.cpp
	   octgrid.h
           my_glut.h)

set(SOURCE_LIGHT main_light.cpp
	   octgrid.cpp
	   octgrid.h)

set(SOURCE_TEST main_test.cpp
	   octgrid.cpp
	   octgrid.h)


add_executable(Octree ${SOURCE})
add_executable(Octree_light ${SOURCE_LIGHT})
add_executable(Octree_test ${SOURCE_TEST})
set(TARGETS Octree Octree_light Octree_test)

foreach(vtarget ${TARGETS})
target_link_libraries(${vtarget} inmost)
#target_link_libraries(${vtarget} lapack blas) 


if(USE_MPI)
  message("linking ${vtarget} with mpi")
  target_link_libraries(${vtarget} ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    message("adding link flags to ${vtarget}")
    set_target_properties(${vtarget} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)

find_package(OpenGL)
find_package(GLUT)

if(OPENGL_FOUND)
	if(GLUT_FOUND)
		include_directories(${OPENGL_INCLUDE_DIR})
		include_directories(${GLUT_INCLUDE_DIR})
		add_definitions(-D__GRAPHICS__)
		target_link_libraries(${vtarget} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
	else(GLUT_FOUND)
		message("GLUT not found for ${vtarget}")
	endif(GLUT_FOUND)
else(OPENGL_FOUND)
	message("OpenGL not found for ${vtarget}")
	
endif(OPENGL_FOUND)

if(USE_PARTITIONER_ZOLTAN)
  message("linking ${vtarget} with Zoltan")
  target_link_libraries(${vtarget} ${ZOLTAN_LIBRARIES})
endif()
if(USE_PARTITIONER_PARMETIS)
  message("linking ${vtarget} with ParMETIS")
  target_link_libraries(${vtarget} ${PARMETIS_LIBRARIES})
endif()

endforeach(vtarget ${TARGETS})

if(NOT GLUT_FOUND)
	set_target_properties(Octree PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)
endif(NOT GLUT_FOUND)
