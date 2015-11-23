#!/bin/bash
#PBS -q regular
#PBS -l nodes=8:ppn=8
#PBS -l walltime=24:00:00
cd ~/MSPP_new/examples/MatSolve/
OPTIONS=petsc_options_parasails.txt
OUTFILE=~/public_html/parasails.txt
for nproc in 64 32 16 8 4 2 1
do
	mpiexec -np $nproc ./main $OPTIONS matrices/A-002.mtx matrices/B-002.rhs | tee -a $OUTFILE
	mpiexec -np $nproc ./main $OPTIONS matrices/A-010.mtx matrices/B-010.rhs | tee -a $OUTFILE
	mpiexec -np $nproc ./main $OPTIONS matrices/A-050.mtx matrices/B-050.rhs | tee -a $OUTFILE
	mpiexec -np $nproc ./main $OPTIONS matrices/2P_320K_3layers_A.mtx matrices/2P_320K_3layers_B.rhs | tee -a $OUTFILE
	mpiexec -np $nproc ./main $OPTIONS matrices/O_320K_3layers_A.mtx matrices/O_320K_3layers_B.rhs | tee -a $OUTFILE
	mpiexec -np $nproc ./main $OPTIONS matrices/2P_GHK_300K_A.mtx matrices/2P_GHK_300K_B.rhs | tee -a $OUTFILE
	mpiexec -np $nproc ./main $OPTIONS matrices/O_GHK_300K_A.mtx matrices/O_GHK_300K_B.rhs | tee -a $OUTFILE
done
#mpiexec gdb ./program -batch -readnow -x commands.gdb
