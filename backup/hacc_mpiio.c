/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Simple MPI-IO program to measure the performance
 * of interleaved pattern vs contiguous pattern.
 * KY
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/stat.h>
#include <mpi.h>
#include <math.h>
#include <stdbool.h>

#define MB          (1024*1024)
#define FILENAME "./data.mpi"
#define NUM_VARS 9

int64_t BUF_SIZE_PER_VAR = 1 * MB;
int64_t NUM_DOUBLES_PER_VAR_PER_RANK;


int main(int argc, char **argv)
{
    int  mpi_size, mpi_rank, i;

    MPI_File fh;
    MPI_Offset  mpi_off = 0;
    MPI_Status  mpi_stat;
    MPI_Info info;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    MPI_Info_create(&info);
    MPI_Info_set(info,"romio_ds_read","disable");
    MPI_Info_set(info,"romio_ds_write","disable");

    BUF_SIZE_PER_VAR = atoi(argv[1]) * MB;
    NUM_DOUBLES_PER_VAR_PER_RANK = BUF_SIZE_PER_VAR / mpi_size / sizeof(double);
    bool collective = (getenv("HACC_CHEN_COLLECTIVE") != NULL);

    // Allocate total amount of data per process to the buffer
    double *writedata = (double*) malloc(BUF_SIZE_PER_VAR/mpi_size*NUM_VARS);
    for (i=0; i < NUM_DOUBLES_PER_VAR_PER_RANK*NUM_VARS; i++) {
        writedata[i] = (double)(mpi_rank*BUF_SIZE_PER_VAR/mpi_size+i);
    }

    MPI_File_open(MPI_COMM_WORLD, FILENAME, MPI_MODE_CREATE|MPI_MODE_WRONLY, info, &fh);

    // write a 2048 header to simulate the HDF5 superblock
    double *garbage = malloc(2048);
    memset(garbage, 100, 2048);
    if(mpi_rank == 0)
        MPI_File_write_at(fh, 0, garbage, 2048, MPI_BYTE, &mpi_stat);

    // Each process has a contiuous write.
    double write_tstart = MPI_Wtime();
    mpi_off = 2048+BUF_SIZE_PER_VAR/mpi_size * mpi_rank; // position of the first variable
    for (i=0; i < NUM_VARS; i++) {
        if(collective) {
            //MPI_File_write_at_all(fh, mpi_off, &writedata[i*NUM_DOUBLES_PER_VAR_PER_RANK], NUM_DOUBLES_PER_VAR_PER_RANK, MPI_DOUBLE, &mpi_stat);
            MPI_File_write_at_all(fh, mpi_off, &writedata[i*NUM_DOUBLES_PER_VAR_PER_RANK], sizeof(double)*NUM_DOUBLES_PER_VAR_PER_RANK, MPI_BYTE, &mpi_stat);
        } else {
            //MPI_File_write_at(fh, mpi_off, &writedata[i*NUM_DOUBLES_PER_VAR_PER_RANK], NUM_DOUBLES_PER_VAR_PER_RANK, MPI_DOUBLE, &mpi_stat);
            MPI_File_write_at(fh, mpi_off, &writedata[i*NUM_DOUBLES_PER_VAR_PER_RANK], sizeof(double)*NUM_DOUBLES_PER_VAR_PER_RANK, MPI_BYTE, &mpi_stat);
        }
        mpi_off += BUF_SIZE_PER_VAR;                // move to the position of next variable
    }
    MPI_File_close(&fh);
    double write_tend = MPI_Wtime();
    free(writedata);

    double min_tstart, max_tend;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&write_tstart, &min_tstart, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&write_tend, &max_tend, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(mpi_rank == 0) {
        double total_time = max_tend - min_tstart;
        printf("Total time: %f, Aggragated Write Bandwidth: %fMB/s\n", total_time, BUF_SIZE_PER_VAR*NUM_VARS/MB/total_time);
    }




    /*
    MPI_File_open(MPI_COMM_WORLD, FILENAME, MPI_MODE_RDONLY, info, &fh);

    // Allocate total amount of data per process to the buffer
    readdata = malloc(buf_size*NUM_VARS/mpi_size*sizeof(double));

    // Each process has a contiuous write.
    for (i=0; i < NUM_VARS; i++) {
        MPI_File_read_at(fh, mpi_off, readdata, buf_size_per_proc*sizeof(double), MPI_BYTE, &mpi_stat);
        mpi_off+=buf_size_per_proc*sizeof(double);
    }


    free(readdata);
    MPI_File_close(&fh);
    */


    MPI_Finalize();

    return 0;
}
