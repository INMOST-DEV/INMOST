CXX=mpicxx
CC=mpicxx


#links for cluster
PETSC_DIR=#/data4t/terekhov/Packages/petsc-3.3-p3
PETSC_ARCH=#arch-linux2-c-opt
INCPATH = #-I/data4t/terekhov/program/include -I${PETSC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include


#OPTFLAGS= -O0 -g -Wall
OPTFLAGS=  -O3

MYLIBS=../../mspp.a
#MYLIBS= ../../msppd.a
CXXFLAGS=$(OPTFLAGS) $(INCPATH)

LDFLAGS=$(OPTFLAGS)
ILU2 = ../../ILU2/lib/libilu-2.3.a 
PETSC = -L${PETSC_DIR}/${PETSC_ARCH}/lib/ 
#-lpetsc -L/usr/X11R6/lib -lX11 -ldmumps -lmumps_common -lmpi_f77 -lscalapack -lpord -lblacs -lparmetis -lmetis -lflapack -lfblas -lmpi -lHYPRE -lsuperlu_dist_2.4
PETSC+=-lpetsc -L/usr/X11R6/lib64 -lX11 -lpthread -lstdc++ -ldmumps -lmumps_common -lscalapack -lpord -lblacs -lparmetis -lmetis -lmpi -lHYPRE -lflapack -lfblas
LDLIBS=$(MYLIBS)  $(PETSC) $(ILU2)


targets=main

all: $(targets)


main: main.o $(MYLIBS)


clean:
	rm -f $(targets) *.o
