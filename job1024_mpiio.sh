#!/bin/bash
#SBATCH -N 32
#SBATCH -n 1024
#SBATCH -t 00:25:00
#SBATCH -p pbatch
#SBATCH --job-name="HACC"

export I_MPI_EXTRA_FILESYSTEM=on
export I_MPI_EXTRA_FILESYSTEM_LIST=lustre

#export I_MPI_OFI_PROVIDER_DUMP=on
#export ROMIO_PRINT_HINTS=1



N=10
VARIABLE_SIZE=8192  # 8GB per variable, 72GB total

: '
echo "===================INDEPENDENT===================="
unset HACC_CHEN_COLLECTIVE
for stripe_size in 36M 72M 128M 512M
do
    echo "stripe size: " $stripe_size
    rm -f ./data.mpi
    lfs setstripe -c 32 -S $stripe_size ./data.mpi

    for i in $(seq 1 $N)
    do
        srun ./hacc_mpiio.out $VARIABLE_SIZE
    done
done


echo "===================COLLECTIVE===================="
export HACC_CHEN_COLLECTIVE=true
for stripe_size in 36M 72M 128M 512M
do
    echo "stripe size: " $stripe_size
    rm -f ./data.mpi
    lfs setstripe -c 32 -S $stripe_size ./data.mpi

    for i in $(seq 1 $N)
    do
        srun ./hacc_mpiio.out $VARIABLE_SIZE
    done
done
'


unset HACC_CHEN_COLLECTIVE

for i in $(seq 1 $N)
do
    rm -f ./data.mpi
    lfs setstripe -c 32 -S 128M ./data.mpi
    srun ./hacc_mpiio.out -c $VARIABLE_SIZE
done

for i in $(seq 1 $N)
do
    rm -f ./data.mpi
    lfs setstripe -c 32 -S 128M ./data.mpi
    srun ./hacc_mpiio.out -i $VARIABLE_SIZE
done



#LD_PRELOAD=/g/g90/wang116/sources/Recorder/lib/librecorder.so srun ./hacc_mpiio.out $VARIABLE_SIZE


#export HACC_CHEN_COLLECTIVE=true
#LD_PRELOAD=/g/g90/wang116/sources/Recorder/lib/librecorder.so srun ./hacc_mpiio.out $VARIABLE_SIZE
#mv logs ./trans/mpiio_collective

