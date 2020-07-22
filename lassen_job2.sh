#!/bin/bash

N=10

nodes=2
tasks_per_node=40
let VAR_SIZE=$nodes*$tasks_per_node*8

echo $nodes, $VAR_SIZE

cd /p/gpfs1/wang116/sources/HACC-IO-HDF5

echo "mpi -i"
for i in $(seq 1 $N)
do
    rm -f ./data.mpi
    lrun -N $nodes -T $tasks_per_node ./hacc_mpiio.out -i $VAR_SIZE
done

echo "hdf5 -m"
for i in $(seq 1 $N)
do
    rm -f ./data.h5
    lrun -N $nodes -T $tasks_per_node ./hacc_hdf5.out -m $VAR_SIZE
done

: '
echo "mpi -c"
for i in $(seq 1 $N)
do
    rm -f ./data.mpi
    lrun -N $nodes -T $tasks_per_node ./hacc_mpiio.out -c $VAR_SIZE
done

echo "hdf5 -i"
for i in $(seq 1 $N)
do
    rm -f ./data.h5
    lrun -N $nodes -T $tasks_per_node ./hacc_hdf5.out -i $VAR_SIZE
done

#export HACC_CHEN_ALIGNMENT=8192
'
