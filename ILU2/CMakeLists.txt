cmake_minimum_required(VERSION 2.8.7)
project(ani3d)
set(SOURCE bcg.f
           matvecCSR.f
           dsort.f
           ilu0.f
           iluoo.f)
#if(NOT WIN32 AND NOT CMAKE_GENERATOR STREQUAL Xcode AND NOT CMAKE_GENERATOR STREQUAL Ninja)
  #--Work-around for CMake issue 0009220
  if(DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER MATCHES "^$")
    set(CMAKE_Fortran_COMPILER CMAKE_Fortran_COMPILER-NOTFOUND)
  endif()
  enable_language(Fortran OPTIONAL)
#endif()
#tutorial for better fortran support:
#http://www.netlib.org/lapack/lawnspdf/lawn270.pdf



if(CMAKE_Fortran_COMPILER_WORKS)
	message("Fortran found, adding ani3d")
	add_library(ani3d STATIC ${SOURCE})
	set(ANI3D_AVAILIBLE 1)
else()
	message("Fortran not found")
	set(ANI3D_AVAILIBLE 0)
endif()
