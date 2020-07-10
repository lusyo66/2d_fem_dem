#!/bin/sh

# #PBS -m abe
# #PBS -M beichuan.yan@colorado.edu
#PBS -j oe

#PBS -N run
#PBS -l nodes=1:ppn=1
#PBS -l walltime=24:00:00

module purge
#module load openmpi-x86_64
module load openmpi-1.6.4-gcc-4.6.4

cd $PBS_O_WORKDIR
mpirun -np 1 ./paraEllip3d_Bagi input.txt

