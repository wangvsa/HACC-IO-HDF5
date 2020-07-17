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
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/stat.h>
#include <stdbool.h>
#include <mpi.h>
#include <math.h>
#include <stdbool.h>
#include "hdf5.h"

#define RANK        1
#define COLL_META   0
#define KB          (1*1024)
#define MB          (1024*1024)



int64_t BUF_SIZE_PER_VAR;                               // Total file size for one variable
int64_t NUM_DOUBLES_PER_VAR;                            // Number of doubles for one variable

int NUM_VARS  = 9;                                      // Number of variables, currently is 9 like Generic IO.
char *filename = "./data.h5";

double read_tstart, read_tend, write_tstart, write_tend;// per process local read/write time


const char * DATASETNAME[] = {
    "id",
    "mask",
    "x",
    "y",
    "z",
    "vx",
    "vy",
    "vz",
    "phi"
};

typedef struct {
    double id;
    double mask;
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;;
    double phi;
} hacc_t;

hid_t get_dcpl_id(int mpi_size) {
    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);

    H5Pset_alloc_time(dcpl_id, H5D_ALLOC_TIME_EARLY);
    H5Pset_fill_time(dcpl_id, H5D_FILL_TIME_NEVER);

    const char* chunk = (getenv("HACC_CHEN_CHUNK_DIM"));
    if(chunk) {
        hsize_t chunk_dims[1] = { atoi(chunk) };
        H5Pset_chunk(dcpl_id, 1, chunk_dims);
    }
    return dcpl_id;
}

hid_t get_fapl_id() {
    hid_t fapl_id = H5Pcreate (H5P_FILE_ACCESS);
    hid_t mpio_fapl_id = H5Pcreate (H5P_FILE_ACCESS);

    const char* align = (getenv("HACC_CHEN_ALIGNMENT"));
    if(align)
        H5Pset_alignment(fapl_id, 0, atoi(align)*KB);

    // H5Pset_dxpl_mpio(fapl_id, H5FD_MPIO_COLLECTIVE);
    if(COLL_META) {
        H5Pset_coll_metadata_write(fapl_id, 1);
        H5Pset_all_coll_metadata_ops(fapl_id, 1 );
    }

    H5Pset_meta_block_size(fapl_id, 8*MB);
    // Use the latest file format
    H5Pset_libver_bounds(fapl_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);

    MPI_Info info;
    MPI_Info_create(&info);
    //MPI_Info_set(info, "romio_cb_read", "enable");
    //MPI_Info_set(info, "romio_cb_write", "enable");
    MPI_Info_set(info, "romio_ds_read", "disable");
    MPI_Info_set(info, "romio_ds_write", "disable");
    MPI_Info_set(info, "romio_lustre_ds_in_coll", "disable");
    //MPI_Info_set(info, "ind_rd_buffer_size", "33554432");
    //MPI_Info_set(info, "cb_buffer_size", "536870912");
    H5Pset_fapl_mpio(mpio_fapl_id, MPI_COMM_WORLD, info);
    H5Pset_fapl_split(fapl_id, "-meta.h5", mpio_fapl_id, "-raw.h5", mpio_fapl_id);
    return fapl_id;
}


hid_t open_hdf5_file(char* filename) {
    hid_t fapl_id = get_fapl_id();
    printf("here 0\n");
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDWR, fapl_id);
    printf("here 1\n");
    H5Pclose(fapl_id);
    return file_id;
}

hid_t create_hdf5_file(char* filename) {
    hid_t file_id;                              /* File ID */
    hid_t fapl_id, fcpl_id;		                /* File access property list */

    fcpl_id = H5Pcreate(H5P_FILE_CREATE);

    // default is H5F_FSPACE_STRATEGY_FSM_AGGR
    // H5F_FSPACE_STRATEGY_PAGE, H5F_FSPACE_STRATEGY_NONE, H5F_FSPACE_STRATEGY_AGGR,
    //H5Pset_file_space_strategy(fcpl_id, H5F_FSPACE_STRATEGY_AGGR, 0, (hsize_t)1);

    fapl_id = get_fapl_id();
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, fcpl_id, fapl_id);

    H5Pclose(fcpl_id);
    H5Pclose(fapl_id);

    return file_id;
}

void hacc_hdf5_with_seperate_dataset(int mpi_rank, int mpi_size) {

    hid_t file_id = create_hdf5_file(filename);

    hid_t dset_id;		        /* Dataset ID */
    hid_t dcpl_id;

    hsize_t file_start[RANK];	/* File dataset selection start coordinates (for hyperslab setting) */
    hsize_t file_count[RANK];	/* File dataset selection count coordinates (for hyperslab setting) */
    hsize_t mem_start[RANK];	/* Memory buffer selection start coordinates (for hyperslab setting) */
    hsize_t mem_count[RANK];	/* Memory buffer selection count coordinates (for hyperslab setting) */

    hid_t file_space_id;	    /* File dataspace ID */
    hid_t mem_space_id;		    /* Memory dataspace ID */
    hsize_t file_dims[RANK];   	/* Dataset dimemsion sizes */
    hsize_t mem_dims[RANK];   	/* Memory buffer dimemsion sizes */

    int i;

    dcpl_id = get_dcpl_id(mpi_size);
    file_dims[0] = BUF_SIZE_PER_VAR / sizeof(double);
    file_space_id = H5Screate_simple(1, file_dims, NULL);
    for (i = 0; i < NUM_VARS; i++) {
        dset_id = H5Dcreate2(file_id, DATASETNAME[i], H5T_NATIVE_DOUBLE, file_space_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        H5Dclose(dset_id);
    }
    H5Sclose(file_space_id);
    H5Pclose(dcpl_id);
    H5Fclose(file_id);

    file_id = open_hdf5_file(filename);

    // 2. Initialize the data for writing
    size_t NUM_DOUBLES_PER_VAR_PER_RANK = NUM_DOUBLES_PER_VAR / mpi_size;
    double *writedata = (double*) malloc(BUF_SIZE_PER_VAR * NUM_VARS / mpi_size);
    for (i = 0; i < NUM_DOUBLES_PER_VAR_PER_RANK * NUM_VARS; i++)
        writedata[i] = (double)(mpi_rank * NUM_DOUBLES_PER_VAR_PER_RANK + i);

    // 3. Write the data
    hid_t dxfer_plist_id = get_dxfer_plist_id();
    mem_dims[0] = NUM_DOUBLES_PER_VAR_PER_RANK * NUM_VARS;   // number of doubles for all varaibles per rank
    file_dims[0] = mem_dims[0] * mpi_size;                   // number of doubles for all variables and all ranks

    mem_space_id = H5Screate_simple(1, mem_dims, NULL);
    file_space_id = H5Screate_simple(1, file_dims, NULL);

    MPI_Barrier(MPI_COMM_WORLD);
    write_tstart = MPI_Wtime();
    for (i = 0; i < NUM_VARS; i++) {

        dset_id = H5Dopen(file_id, DATASETNAME[i], H5P_DEFAULT);

        // Select column of elements in the file dataset
        file_start[0] = mpi_rank * NUM_DOUBLES_PER_VAR_PER_RANK;
        file_count[0] = NUM_DOUBLES_PER_VAR_PER_RANK;
        H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

        // Select all elements in the memory buffer
        mem_start[0] = i * NUM_DOUBLES_PER_VAR_PER_RANK;
        mem_count[0] = NUM_DOUBLES_PER_VAR_PER_RANK;
        H5Sselect_hyperslab(mem_space_id, H5S_SELECT_SET, mem_start, NULL, mem_count, NULL);

        // Write data independently
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, mem_space_id, file_space_id, dxfer_plist_id, writedata);

        H5Dclose(dset_id);
    }
    write_tend = MPI_Wtime();
    H5Fflush(file_id, H5F_SCOPE_GLOBAL);

    MPI_Barrier(MPI_COMM_WORLD);
    free(writedata);
}



int main(int argc, char **argv)
{

    int mpi_size, mpi_rank;	    /* MPI variables */

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

}

