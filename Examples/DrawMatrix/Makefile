CXX=mpicxx
CC=mpicxx

PETSC_DIR=#/data4t/terekhov/Packages/petsc-3.3-p3
PETSC_ARCH=#arch-linux2-c-opt

INCPATH = -I../../ -I${PETSC_DIR}/${PETSC_ARCH}/include -I${PETSC_DIR}/include

OPENGL = -L/usr/X11R6/lib -lX11 -lXi -lXmu -lGL -lglut -lGLU
#OPENGL = -framework GLUT -framework OpenGL
#OPENGL = -lglut -lgl -lglu32 -lopengl32
#OPENGL = -L/lib -lglu32 -lopengl32 

PETSC = -L${PETSC_DIR}/${PETSC_ARCH}/lib/ -lpetsc -L/usr/X11R6/lib -lX11 -lparmetis -lmetis -lmpi_f77 -lflapack -lfblas -lgfortran
ILU2 = ../../ILU2/lib/libilu-2.3.a -lgfortran

#OPTFLAGS= -O0 -g -Wall
OPTFLAGS=  -O3

MYLIBS=../../mspp.a
#MYLIBS= ../../msppd.a
CXXFLAGS=$(OPTFLAGS) $(INCPATH)

LDFLAGS=$(OPTFLAGS) 
LDLIBS=$(MYLIBS) $(OPENGL) $(PETSC) $(ILU2)


targets=main

all: $(targets)


main: main.o


clean:
	rm -f $(targets) *.o
