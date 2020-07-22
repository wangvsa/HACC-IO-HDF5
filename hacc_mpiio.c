#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <mpi.h>
#include <math.h>
#include <stdbool.h>

#define KB          (1024L)
#define MB          (1024L*1024L)
#define FILENAME "./data.mpi"
#define NUM_VARS 9

size_t BUF_SIZE_PER_VAR = 1 * MB;
size_t NUM_DOUBLES_PER_VAR_PER_RANK;

void print_io_stat(int rank, int mpi_size) {
    pid_t pid = getpid();

    char filename[256];
    sprintf(filename, "/proc/%d/io", pid);

    char data[1024];
    char rbuf[mpi_size*1024];
    FILE* f = fopen(filename, "r");
    fread(data, sizeof(char), 1024, f);
    fclose(f);

    MPI_Gather(data, 1024, MPI_BYTE, rbuf, 1024, MPI_BYTE, 0, MPI_COMM_WORLD);

    if(rank == 0) {
        f = fopen("./cache_stat.txt", "a");
        int i;
        for(i = 0; i < mpi_size; i++) {
            int off = i * 1024;
            fprintf(f, "\nFrom Rank %d:\n%s", i, rbuf+off);
        }
        fclose(f);
    }
}


void write_contiguous(size_t start_offset, int mpi_rank, int mpi_size)
{
    int i;
    MPI_File fh;
    MPI_Offset  mpi_off = 0;
    MPI_Status  mpi_stat;
    MPI_Info info;

    MPI_Info_create(&info);
    MPI_Info_set(info,"romio_ds_read","disable");
    MPI_Info_set(info,"romio_ds_write","disable");

    NUM_DOUBLES_PER_VAR_PER_RANK = BUF_SIZE_PER_VAR / mpi_size / sizeof(double);
    bool collective = (getenv("HACC_CHEN_COLLECTIVE") != NULL);

    // Allocate total amount of data per process to the buffer
    double *writedata = (double*) malloc(BUF_SIZE_PER_VAR/mpi_size*NUM_VARS);
    for (i=0; i < NUM_DOUBLES_PER_VAR_PER_RANK*NUM_VARS; i++) {
        writedata[i] = (double)(mpi_rank * NUM_DOUBLES_PER_VAR_PER_RANK+ i % NUM_DOUBLES_PER_VAR_PER_RANK);
    }

    MPI_File_open(MPI_COMM_WORLD, FILENAME, MPI_MODE_CREATE|MPI_MODE_WRONLY, info, &fh);

    // Each process has a contiuous write.
    MPI_Barrier(MPI_COMM_WORLD);
    double write_tstart = MPI_Wtime();
    mpi_off = start_offset+BUF_SIZE_PER_VAR/mpi_size * mpi_rank; // position of the first variable
    for (i=0; i < NUM_VARS; i++) {
        if(collective) {
            MPI_File_write_at_all(fh, mpi_off, &writedata[i*NUM_DOUBLES_PER_VAR_PER_RANK], NUM_DOUBLES_PER_VAR_PER_RANK, MPI_DOUBLE, &mpi_stat);
        } else {
            MPI_File_write_at(fh, mpi_off, &writedata[i*NUM_DOUBLES_PER_VAR_PER_RANK], NUM_DOUBLES_PER_VAR_PER_RANK, MPI_DOUBLE, &mpi_stat);
        }
        mpi_off += BUF_SIZE_PER_VAR;                // move to the position of next variable
    }
    // MPI_File_sync(fh);
    double write_tend = MPI_Wtime();
    MPI_File_close(&fh);
    //print_io_stat(mpi_rank, mpi_size);
    free(writedata);

    double min_tstart, max_tend;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&write_tstart, &min_tstart, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&write_tend, &max_tend, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(mpi_rank == 0) {
        double total_time = max_tend - min_tstart;
        printf("Total time: %f, Write Bandwidth: %fMB/s\n", total_time, BUF_SIZE_PER_VAR*NUM_VARS/MB/total_time);
    }
}

/**
 * From the view of variables, each write is interleaved
 * From the view of a single process, 9 writes are contiguous.
 *
 * This has the same pattern as HDF5+Multi
 */
void write_interleaved(size_t start_offset, int mpi_rank, int mpi_size)
{
    int i;
    MPI_File fh;
    MPI_Offset  mpi_off = 0;
    MPI_Status  mpi_stat;
    MPI_Info info;

    MPI_Info_create(&info);
    MPI_Info_set(info,"romio_ds_read","disable");
    MPI_Info_set(info,"romio_ds_write","disable");
    MPI_Info_set(info, "romio_cb_read", "disable");
    MPI_Info_set(info, "romio_cb_write", "disable");


    NUM_DOUBLES_PER_VAR_PER_RANK = BUF_SIZE_PER_VAR / mpi_size / sizeof(double);
    bool collective = (getenv("HACC_CHEN_COLLECTIVE") != NULL);

    // Allocate total amount of data per process to the buffer
    double *writedata = (double*) malloc(BUF_SIZE_PER_VAR/mpi_size*NUM_VARS);
    for (i=0; i < NUM_DOUBLES_PER_VAR_PER_RANK*NUM_VARS; i++) {
        writedata[i] = (double)(mpi_rank * NUM_DOUBLES_PER_VAR_PER_RANK+ i % NUM_DOUBLES_PER_VAR_PER_RANK);
    }

    MPI_File_open(MPI_COMM_WORLD, FILENAME, MPI_MODE_CREATE|MPI_MODE_WRONLY, info, &fh);

    // Each process has a contiuous write.
    MPI_Barrier(MPI_COMM_WORLD);
    double write_tstart = MPI_Wtime();
    mpi_off = start_offset+BUF_SIZE_PER_VAR/mpi_size*NUM_VARS * mpi_rank; // position of the first variable
    for (i=0; i < NUM_VARS; i++) {
        if(collective) {
            MPI_File_write_at_all(fh, mpi_off, &writedata[i*NUM_DOUBLES_PER_VAR_PER_RANK], NUM_DOUBLES_PER_VAR_PER_RANK, MPI_DOUBLE, &mpi_stat);
        } else {
            MPI_File_write_at(fh, mpi_off, &writedata[i*NUM_DOUBLES_PER_VAR_PER_RANK], NUM_DOUBLES_PER_VAR_PER_RANK, MPI_DOUBLE, &mpi_stat);
        }
        mpi_off += BUF_SIZE_PER_VAR/mpi_size;                // move to the position of next variable
    }
    // MPI_File_sync(fh);
    double write_tend = MPI_Wtime();
    MPI_File_close(&fh);
    //print_io_stat(mpi_rank, mpi_size);
    free(writedata);

    double min_tstart, max_tend;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&write_tstart, &min_tstart, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&write_tend, &max_tend, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(mpi_rank == 0) {
        double total_time = max_tend - min_tstart;
        printf("Total time: %f, Write Bandwidth: %fMB/s\n", total_time, BUF_SIZE_PER_VAR*NUM_VARS/MB/total_time);
    }
}

int main(int argc, char **argv)
{

    int mpi_size, mpi_rank;
    int i;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    int interleaved = 0;
    if (strcmp(argv[1], "-i") == 0) interleaved = 1;

    BUF_SIZE_PER_VAR = atoi(argv[2]) * MB;

    size_t offsets[] = {0, 32*MB};
    //if(mpi_rank == 0) printf("start offset: %ldKB\n", start_offset/KB);
    for(i = 0; i < 2; i++) {
        size_t start_offset = offsets[i];
        if(mpi_rank == 0)
            printf("offset: %ld\n", start_offset);

            if(mpi_rank == 0)
                remove(FILENAME);
            if(interleaved)
                write_interleaved(start_offset, mpi_rank, mpi_size);
            else
                write_contiguous(start_offset, mpi_rank, mpi_size);
    }

    MPI_Finalize();
    return 0;
}
