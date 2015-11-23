#ifndef INMOST_OPTIONS_CMAKE_INCLUDED
#define INMOST_OPTIONS_CMAKE_INCLUDED

#cmakedefine USE_OMP

#cmakedefine USE_MESH

#cmakedefine USE_AUTODIFF
#cmakedefine USE_AUTODIFF_EXPRESSION_TEMPLATES
#cmakedefine USE_AUTODIFF_OPENCL
#cmakedefine USE_AUTODIFF_ASMJIT

#cmakedefine USE_PARTITIONER
#cmakedefine USE_PARTITIONER_ZOLTAN
#cmakedefine USE_PARTITIONER_PARMETIS

#cmakedefine USE_SOLVER
#cmakedefine USE_SOLVER_MONDRIAAN
#cmakedefine USE_SOLVER_METIS
#cmakedefine USE_SOLVER_PETSC
#cmakedefine USE_SOLVER_TRILINOS
#cmakedefine USE_SOLVER_ANI

#cmakedefine USE_NONLINEAR
#cmakedefine USE_NONLINEAR_TRILINOS
#cmakedefine USE_NONLINEAR_PETSC
#cmakedefine USE_NONLINEAR_SUNDIALS

#cmakedefine USE_MPI //include mpi for mpi functions
#cmakedefine USE_MPI_P2P //use (probably) more effective point to point algorithms
#cmakedefine USE_MPI_FILE //use functionality for parallel files
#cmakedefine USE_MPI2 //use mpi-2 extensions


#endif //INMOST_OPTIONS_CMAKE_INCLUDED
