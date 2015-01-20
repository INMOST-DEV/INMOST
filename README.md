INMOST
======

A toolkit for distributed mathematical modeling

## New address
INMOST was moved to the new address: https://github.com/INMOST-DEV/INMOST

Please update your bookmarks, if you used `git clone` before then refer to https://github.com/INMOST-DEV/INMOST/issues/7#issuecomment-70378751

## Install Guide
Refer to https://github.com/INMOST-DEV/INMOST/wiki/0100-Compilation

## Examples

Several representative examples are provided in source archive.
Here we will try three parallel steps: grid generation, FVM discretization and linear matrix solution.
Each example may be executed in serial or parallel ways.

### Parallel Grid Generation

This example creates simple cubic or prismatic mesh. You can use ParaView to view the meshes.
```
cd examples/GridGen
mpirun -np 4 ./GridGen 4 32 32 32
```
Generator parameters are: `ng nx ny nz`
where `ng=3` stands for Prismatic generator and
`ng=4` for Cubic one,
while `nx ny nz` are the mesh dimensions.

File grid.pvtk (as well as grid_X.vtk with X=0,1,2,3) will appear in the current directory.
Run
`paraview --data=grid.pvtk`
and try the following tags in objects to display:
- `P_OWNER_PROCESSOR` – partitioning to processors
- `GLOBAL_ID` – global cell ID number

### Parallel Finite Volume Discretization

This example uses simple two-point FVM scheme to solve Laplace's equation in unit cube domain.
```
cd ../FVDiscr
mpirun -np 4 ./FVDiscr ../GridGen/grid.pvtk A.mtx b.rhs
```
Files result.pvtk (as well as result_X.vtk with X=0,1,2,3) and A.mtx b.rhs will appear in the current directory.
Run
`paraview --data=result.pvtk`
and try the following tags in objects to display:
- `Solution` – the solution to the problem
- `K` – tensor K (constant equal to 1 in this example)

### Solve the Matrix stored in mtx format

This example solves the linear system using different solvers.
```
cd ../MatSolve
mpirun -np 4 ./MatSolve 0 ../FVDiscr/A.mtx ../FVDiscr/b.rhs
```
Solution time and the true residual will output to the screen.
The first parameter selects the solver:
- `0` – INMOST ILU2
- `1` – INMOST MLILUC
- `2` – ANI3D
- `3` – PETSc

### Other Examples

There are other examples:
```
DrawGrid
DrawMatrix
OctreeCutcell
OldDrawGrid
Solver
```
which can be used without warranty.

