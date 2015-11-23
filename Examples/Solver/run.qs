#!/bin/bash
#PBS -q sixcore
#PBS -l nodes=6:ppn=12
#PBS -l walltime=24:00:00
cd ~/MSPP_new/examples/Solver/
#mpiexec -np 30 ./main ../Grids/ghk_300k.vtk -ksp_type bcgs -pc_type asm -pc_asm_type basic -pc_asm_overlap 1 -sub_pc_type ilu -sub_pc_factor_levels 1
#cp grid.pmf ../Grid/ghk_300k_30rcm.pmf
#mpiexec -np 18 ./main ../Grids/ghk_300k.vtk -ksp_type bcgs -pc_type asm -pc_asm_type basic -pc_asm_overlap 1 -sub_pc_type ilu -sub_pc_factor_levels 1
#cp grid.pmf ../Grid/ghk_300k_18rcm.pmf

mpiexec -np 72 ./main ../Grids/ghk_300k_72rcm.pmf -ksp_type bcgs -pc_type asm -pc_asm_type basic -pc_asm_overlap 1 -sub_pc_type ilu -sub_pc_factor_levels 1
mpiexec -np 60 ./main ../Grids/ghk_300k_60rcm.pmf -ksp_type bcgs -pc_type asm -pc_asm_type basic -pc_asm_overlap 1 -sub_pc_type ilu -sub_pc_factor_levels 1
mpiexec -np 48 ./main ../Grids/ghk_300k_48rcm.pmf -ksp_type bcgs -pc_type asm -pc_asm_type basic -pc_asm_overlap 1 -sub_pc_type ilu -sub_pc_factor_levels 1
mpiexec -np 36 ./main ../Grids/ghk_300k_36rcm.pmf -ksp_type bcgs -pc_type asm -pc_asm_type basic -pc_asm_overlap 1 -sub_pc_type ilu -sub_pc_factor_levels 1
mpiexec -np 24 ./main ../Grids/ghk_300k_24rcm.pmf -ksp_type bcgs -pc_type asm -pc_asm_type basic -pc_asm_overlap 1 -sub_pc_type ilu -sub_pc_factor_levels 1
mpiexec -np 12 ./main ../Grids/ghk_300k_12rcm.pmf -ksp_type bcgs -pc_type asm -pc_asm_type basic -pc_asm_overlap 1 -sub_pc_type ilu -sub_pc_factor_levels 1


mpiexec -np 72 ./main ../Grids/ghk_300k_36rcm.pmf -ksp_type bcgs -pc_type asm -pc_asm_type basic -pc_asm_overlap 1 -sub_pc_type ilu -sub_pc_factor_levels 1
mpiexec -np 60 ./main ../Grids/ghk_300k_30rcm.pmf -ksp_type bcgs -pc_type asm -pc_asm_type basic -pc_asm_overlap 1 -sub_pc_type ilu -sub_pc_factor_levels 1
mpiexec -np 48 ./main ../Grids/ghk_300k_24rcm.pmf -ksp_type bcgs -pc_type asm -pc_asm_type basic -pc_asm_overlap 1 -sub_pc_type ilu -sub_pc_factor_levels 1
mpiexec -np 36 ./main ../Grids/ghk_300k_18rcm.pmf -ksp_type bcgs -pc_type asm -pc_asm_type basic -pc_asm_overlap 1 -sub_pc_type ilu -sub_pc_factor_levels 1
mpiexec -np 24 ./main ../Grids/ghk_300k_12rcm.pmf -ksp_type bcgs -pc_type asm -pc_asm_type basic -pc_asm_overlap 1 -sub_pc_type ilu -sub_pc_factor_levels 1
mpiexec -np 12 ./main ../Grids/ghk_300k_6rcm.pmf -ksp_type bcgs -pc_type asm -pc_asm_type basic -pc_asm_overlap 1 -sub_pc_type ilu -sub_pc_factor_levels 1

#mpiexec gdb ./program -batch -readnow -x commands.gdb
