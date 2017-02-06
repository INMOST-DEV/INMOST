project(Octree)
set(SOURCE main.cpp
           rotate.cpp
           rotate.h
	   octgrid.cpp
	   octgrid.h
           my_glut.h)

add_executable(Octree ${SOURCE})

target_link_libraries(Octree inmost lapack blas) 


if(USE_MPI)
  message("linking with mpi")
  target_link_libraries(Octree ${MPI_LIBRARIES}) 
  if(MPI_LINK_FLAGS)
    message("adding link flags")
    set_target_properties(Octree PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif() 
endif(USE_MPI)

find_package(OpenGL)
find_package(GLUT)

if(OPENGL_FOUND)
	if(GLUT_FOUND)
		include_directories(${OPENGL_INCLUDE_DIR})
		include_directories(${GLUT_INCLUDE_DIR})
		add_definitions(-D__GRAPHICS__ -g)
		target_link_libraries(Octree ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
	else(GLUT_FOUND)
		message("GLUT not found")
	endif(GLUT_FOUND)
else(OPENGL_FOUND)
	message("OpenGL not found")
endif(OPENGL_FOUND)

if(USE_PARTITIONER_ZOLTAN)
  message("linking Octree with Zoltan")
  target_link_libraries(Octree ${ZOLTAN_LIBRARIES})
endif()
if(USE_PARTITIONER_PARMETIS)
  message("linking Octree with ParMETIS")
  target_link_libraries(Octree ${PARMETIS_LIBRARIES})
endif()

