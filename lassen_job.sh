#!/bin/bash

### LSF syntax

#BSUB -nnodes 32                  #number of nodes
#BSUB -W 30                       #walltime in minutes
#BSUB -q pbatch                   #queue to use

N=10
VAR_SIZE=8192

cd /p/gpfs1/wang116/sources/HACC-IO-HDF5

: '
echo "mpi -i"
for i in $(seq 1 $N)
do
    rm -f ./data.mpi
    lrun -n 1024 ./hacc_mpiio.out -i 8192
done

echo "mpi -c"
for i in $(seq 1 $N)
do
    rm -f ./data.mpi
    lrun -n 1024 ./hacc_mpiio.out -c 8192
done
'

echo "hdf5 -i"
for i in $(seq 1 $N)
do
    rm -f ./data.h5
    lrun -n 1024 ./hacc_hdf5.out -i 8192
done

echo "hdf5 -c"
for i in $(seq 1 $N)
do
    rm -f ./data.h5
    lrun -n 1024 ./hacc_hdf5.out -c 8192
done

echo "hdf5 -m"
for i in $(seq 1 $N)
do
    rm -f ./data.h5
    lrun -n 1024 ./hacc_hdf5.out -m 8192
done
