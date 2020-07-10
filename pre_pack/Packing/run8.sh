#!/bin/sh

# #PBS -m abe
# #PBS -M beichuan.yan@colorado.edu
#PBS -j oe

#PBS -N run
#PBS -l nodes=1:ppn=8
#PBS -l walltime=100:00:00

module purge
#module load openmpi-x86_64
module load openmpi-4.0.2-gcc-9.2.0

export LD_LIBRARY_PATH=/usr/local/boost-1.70.0-openmpi-4.0.2-gcc-9.2.0/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/qhull-2015.2-gcc-9.2.0/lib:$LD_LIBRARY_PATH


#export LD_LIBRARY_PATH=/usr/local/qhull-2015.2/lib:$LD_LIBRARY_PATH

cd $PBS_O_WORKDIR
mpirun -np 8 ./paraEllip3d input.txt

