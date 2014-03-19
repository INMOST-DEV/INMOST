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
#cmakedefine USE_SOLVER_PETSC
#cmakedefine USE_SOLVER_ANI


#cmakedefine USE_MPI //include mpi for mpi functions
#cmakedefine USE_MPI2 //use (probably) more effective mpi-2 algorithms


#endif //INMOST_OPTIONS_CMAKE_INCLUDED
