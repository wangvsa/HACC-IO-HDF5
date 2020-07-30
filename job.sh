#!/bin/bash
#SBATCH -N 4
#SBATCH -n 128
#SBATCH -t 00:10:00
#SBATCH -p pdebug
#SBATCH --job-name="HACC"

export I_MPI_EXTRA_FILESYSTEM=on
export I_MPI_EXTRA_FILESYSTEM_LIST=lustre

#export I_MPI_DEBUG=2
#export I_MPI_STATS=20
#export I_MPI_OFI_PROVIDER_DUMP=on
#export ROMIO_PRINT_HINTS=1


N=10
VARIABLE_SIZE=1024

rm -f ./data.h5*
lfs setstripe -c 4 -S 128M ./data.mpi
srun ./hacc_mpiio.out -i $VARIABLE_SIZE
#LD_PRELOAD=/g/g90/wang116/sources/Recorder/lib/librecorder.so srun ./hacc_hdf5.out -i $VARIABLE_SIZE


: '
# Try to find the best stripe size
for stripe_size in 1M 2M 4M 8M 16M 32M 64M 128M
do
    echo "stripe size: " $stripe_size
    rm -f ./data.h5
    lfs setstripe -c 4 -S $stripe_size ./data.h5

    echo "Write with compound datatype"
    for i in $(seq 1 $N)
    do
        srun ./hacc_hdf5.out -c $VARIABLE_SIZE
    done
    echo "============================="

    echo "Write with individual dataset"
    for i in $(seq 1 $N)
    do
        srun ./hacc_hdf5.out -i $VARIABLE_SIZE
    done
    echo "============================="

done
'

# Try to find the best stripe count
: '
for stripe_count in 1 2 4 8 16 32 64 96 128
do
    echo "stripe count: " $stripe_count
    rm -f ./data.h5
    lfs setstripe -S 4M -c $stripe_count ./data.h5

    echo "Write with compound datatype"
    for i in $(seq 1 $N)
    do
        srun ./hacc_hdf5.out -c $VARIABLE_SIZE
    done
    echo "============================="

    echo "Write with individual dataset"
    for i in $(seq 1 $N)
    do
        srun ./hacc_hdf5.out -i $VARIABLE_SIZE
    done
    echo "============================="

done
'

: '
# Find the best alignment size
rm -f ./data.h5
lfs setstripe -c 4 -S 32M ./data.h5
for alignment in # 2 4 32 128 4096
do
    #export HACC_CHEN_ALIGNMENT=$alignment
    echo "Write with individual dataset"
    for i in $(seq 1 $N)
    do
        srun ./hacc_hdf5.out -i $VARIABLE_SIZE
    done

    echo "Write with compound datatype" $HACC_CHEN_ALIGNMENT
    for i in $(seq 1 $N)
    do
        srun ./hacc_hdf5.out -c $VARIABLE_SIZE
    done
done
'

# Find the best chunk size
# each dataset has 100*1024*1024/8 = 13107200 doubles
: '
rm -f ./data.h5
lfs setstripe -c 4 -S 32M ./data.h5
for chunk in 819200, 1638400, 3276800, 6553600, 13107200
do
    export HACC_CHEN_CHUNK_DIM=$chunk
    echo "Write with compound datatype" $HACC_CHEN_CHUNK_DIM
    for i in $(seq 1 $N)
    do
        srun ./hacc_hdf5.out -c $VARIABLE_SIZE
    done

    echo "Write with individual dataset"
    for i in $(seq 1 $N)
    do
        srun ./hacc_hdf5.out -i $VARIABLE_SIZE
    done
done
'
