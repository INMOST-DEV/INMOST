#!/bin/bash
#PBS -N MatSol
#PBS -q x12core
#PBS -l nodes=1:ppn=12
#PBS -l walltime=0:30:00
cd $PBS_O_WORKDIR
dir=`basename $PBS_O_WORKDIR`
rm -f res

#mat="z32k.mtx z32k.rhs - database.txt"
#mat="z32k.mtx - 1 database.txt"
mat="n142-2.mtx n142-2.rhs - database.txt"

for p in 4 ; do ### 1 2 3 4 5 6 7 8 9 10 11 12  12 24 36 72 144
#for s in 8 0 ; do ### 8-PE 10-FC 11-K3 0-II
for s in 8 10 11 0 ; do ### 8-PE 10-FC 11-K3 0-II

if [ $p -le 12 ] ; then npernode=$p ; else npernode=12 ; fi

echo ::: dir=$dir p=$p npernode=$npernode s=$s mat=$mat ::: >> res
mpiexec -np $p -npernode $npernode ./MatSolve $s $mat >> res

done; done
