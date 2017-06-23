cmake_minimum_required(VERSION 2.8.7)
project(ani3d)
set(SOURCE inmost_ani_fem.cpp
           inmost_ani_mesh.cpp
           utils.cpp)
  if(DEFINED CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER MATCHES "^$")
    set(CMAKE_Fortran_COMPILER CMAKE_Fortran_COMPILER-NOTFOUND)
  endif()
  enable_language(Fortran OPTIONAL)
#tutorial for better fortran support:
#http://www.netlib.org/lapack/lawnspdf/lawn270.pdf

add_library(ani_inmost STATIC ${SOURCE})


