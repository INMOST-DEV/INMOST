############################################################
# Fortran and C compilers for SEQuential and PARallel codes.
############################################################
F77     = gfortran -m64 #compilers
CC 	= gcc
CXX     = g++
VIEWER  = gmv      #visualization of GMV files

FLINKER = $(F77)   #linkers
CLINKER = $(CC)
LD      = $(CXX)  

LIBSYS  = -lgfortran        #for linking C and Fortran


############################################################
FFLAGS  = -O5 #-fopenmp  #-Wall
CFLAGS  = -O5   
LDFLAGS = 


############################################################

ANIHOME = $(SYSHOME)

ANIILU  = $(ANIHOME)/src/aniILU 

ANILIB 	= $(ANIHOME)/lib
ANIINC  = $(ANIHOME)/include

ANIBLAS   = $(ANIHOME)/blas
ANILAPACK = $(ANIHOME)/src/lapack


############################################################
LIBILU  = $(ANILIB)/libilu-2.3.a 

#LIBBLAS   = -lblas
#LIBLAPACK = -llapack
LIBBLAS   = $(ANILIB)/libblas-3.0.a
LIBLAPACK = $(ANILIB)/liblapack-3.0.a


############################################################
INCLUDE = 



