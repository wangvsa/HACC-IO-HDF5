#!/bin/bash
#SBATCH -N 32
#SBATCH -n 1024
#SBATCH -t 00:15:00
#SBATCH -p pbatch
#SBATCH --job-name="HACC"

export I_MPI_EXTRA_FILESYSTEM=on
export I_MPI_EXTRA_FILESYSTEM_LIST=lustre

#export I_MPI_OFI_PROVIDER_DUMP=on
#export ROMIO_PRINT_HINTS=1



N=10
VARIABLE_SIZE=8192 # 8GB per variable, 72GB total


echo "===================INDEPENDENT===================="
unset HACC_CHEN_COLLECTIVE

for stripe_size in 128M
do
    echo "72GB, Stripe Count 32, Stripe Size" $stripe_size
    #export HACC_CHEN_ALIGNMENT=$alignment
    echo "-i"
    for i in $(seq 1 $N)
    do
        rm -f ./data.h5
        lfs setstripe -c 32 -S $stripe_size ./data.h5
        srun ./hacc_hdf5.out -i $VARIABLE_SIZE
    done

    echo "-c"
    for i in $(seq 1 $N)
    do
        rm -f ./data.h5
        lfs setstripe -c 32 -S $stripe_size ./data.h5
        srun ./hacc_hdf5.out -c $VARIABLE_SIZE
    done

    echo "-m"
    for i in $(seq 1 $N)
    do
        rm -f ./data.h5
        lfs setstripe -c 32 -S $stripe_size ./data.h5
        srun ./hacc_hdf5.out -m $VARIABLE_SIZE
    done
done

: '
echo "===================COLLECTIVE===================="
export HACC_CHEN_COLLECTIVE=true
echo "-c"
for i in $(seq 1 $N)
do
    srun ./hacc_hdf5.out -c $VARIABLE_SIZE
done

echo "-i"
for i in $(seq 1 $N)
do
    srun ./hacc_hdf5.out -i $VARIABLE_SIZE
done

echo "-m"
for i in $(seq 1 $N)
do
    srun ./hacc_hdf5.out -m $VARIABLE_SIZE
done
'
