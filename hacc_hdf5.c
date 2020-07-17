#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdbool.h>
#include <mpi.h>
#include <math.h>
#include <stdbool.h>
#include "hdf5.h"

#define RANK        1
#define COLL_META   0
#define KB          (1*1024L)
#define MB          (1024L*1024L)



size_t BUF_SIZE_PER_VAR;                               // Total file size for one variable
size_t NUM_DOUBLES_PER_VAR;                            // Number of doubles for one variable

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

    //H5Pset_meta_block_size(fapl_id, 8*MB);
    // Use the latest file format
    //H5Pset_libver_bounds(fapl_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);

    MPI_Info info;
    MPI_Info_create(&info);
    //MPI_Info_set(info, "romio_cb_read", "enable");
    //MPI_Info_set(info, "romio_cb_write", "enable");
    MPI_Info_set(info, "romio_ds_read", "disable");
    MPI_Info_set(info, "romio_ds_write", "disable");
    //MPI_Info_set(info, "romio_lustre_ds_in_coll", "disable");
    //MPI_Info_set(info, "ind_rd_buffer_size", "33554432");
    //MPI_Info_set(info, "cb_buffer_size", "536870912");
    H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, info);
    //H5Pset_fapl_split(fapl_id, "-meta.h5", mpio_fapl_id, "-raw.h5", mpio_fapl_id);
    return fapl_id;
}

hid_t get_dxfer_plist_id() {
    const char* collective = (getenv("HACC_CHEN_COLLECTIVE"));
    if(collective) {
        hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        return plist_id;
    } else {
        return H5P_DEFAULT;
    }
}

hid_t open_hdf5_file(char* filename) {
    hid_t fapl_id = get_fapl_id();
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDWR, fapl_id);
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

void hacc_hdf5_with_compound_type(int mpi_rank, int mpi_size) {

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

    // 1. Create a compound datatype to store all 9 variables
    hacc_t *Hdata;
    hid_t Hmemtype;

    Hmemtype = H5Tcreate (H5T_COMPOUND, sizeof (hacc_t));
    H5Tinsert (Hmemtype, "id", HOFFSET (hacc_t, id), H5T_NATIVE_DOUBLE);
    H5Tinsert (Hmemtype, "mask", HOFFSET (hacc_t, mask), H5T_NATIVE_DOUBLE);
    H5Tinsert (Hmemtype, "x", HOFFSET (hacc_t, x), H5T_NATIVE_DOUBLE);
    H5Tinsert (Hmemtype, "y", HOFFSET (hacc_t, y), H5T_NATIVE_DOUBLE);
    H5Tinsert (Hmemtype, "z", HOFFSET (hacc_t, z), H5T_NATIVE_DOUBLE);
    H5Tinsert (Hmemtype, "vx", HOFFSET (hacc_t, vx), H5T_NATIVE_DOUBLE);
    H5Tinsert (Hmemtype, "vy", HOFFSET (hacc_t, vy), H5T_NATIVE_DOUBLE);
    H5Tinsert (Hmemtype, "vz", HOFFSET (hacc_t, vz), H5T_NATIVE_DOUBLE);
    H5Tinsert (Hmemtype, "phi", HOFFSET (hacc_t, phi), H5T_NATIVE_DOUBLE);

    // 2. Initialize the data for writing
    size_t BUF_SIZE_PER_RANK = BUF_SIZE_PER_VAR * NUM_VARS / mpi_size;
    size_t NUM_HDATA_PER_RANK = BUF_SIZE_PER_RANK / sizeof(hacc_t);

    Hdata = (hacc_t *) malloc(BUF_SIZE_PER_RANK);
    size_t i;
    for (i=0; i < NUM_HDATA_PER_RANK; i++) {
        Hdata[i].id = (double)(mpi_rank*NUM_HDATA_PER_RANK+ i);
        Hdata[i].mask =  Hdata[i].id;
        Hdata[i].x =  Hdata[i].id;
        Hdata[i].y =  Hdata[i].id;
        Hdata[i].z =  Hdata[i].id;
        Hdata[i].vx =  Hdata[i].id;
        Hdata[i].vy =  Hdata[i].id;
        Hdata[i].vz =  Hdata[i].id;
        Hdata[i].phi =  Hdata[i].id;
    }

    // 3. Create the dataset
    mem_dims[0] = NUM_HDATA_PER_RANK;
    file_dims[0] = mem_dims[0]*mpi_size;
    file_space_id = H5Screate_simple(1, file_dims, NULL);
    dcpl_id = get_dcpl_id(mpi_size);
    dset_id = H5Dcreate2(file_id, "ALLVAR", Hmemtype, file_space_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    H5Sclose(file_space_id);
    H5Pclose(dcpl_id);
    H5Dclose(dset_id);

    // 3.x prepare for read/write
    file_space_id = H5Screate_simple(1, file_dims, NULL);
    mem_space_id = H5Screate_simple(1, mem_dims, NULL);
    hid_t dxfer_plist_id = get_dxfer_plist_id();

    // 4. Write the data
    MPI_Barrier(MPI_COMM_WORLD);
    write_tstart = MPI_Wtime();
    dset_id = H5Dopen(file_id, "ALLVAR", H5P_DEFAULT);

    // Select column of elements in the file dataset
    file_start[0] = mpi_rank*mem_dims[0];
    file_count[0] = mem_dims[0];
    H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

    // Select all elements in the memory buffer
    mem_start[0] = 0;
    mem_count[0] = mem_dims[0];
    H5Sselect_hyperslab(mem_space_id, H5S_SELECT_SET, mem_start, NULL, mem_count, NULL);

    // Write data independently
    H5Dwrite(dset_id, Hmemtype, mem_space_id, file_space_id, dxfer_plist_id, Hdata);
    H5Dclose(dset_id);
    write_tend = MPI_Wtime();
    H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    free(Hdata);

    // 5. Read data
    hacc_t* readdata = (hacc_t*) malloc(BUF_SIZE_PER_RANK);
    read_tstart = MPI_Wtime();
    dset_id = H5Dopen(file_id, "ALLVAR", H5P_DEFAULT);

    // Select column of elements in the file dataset
    file_start[0] = mpi_rank*mem_dims[0];
    file_count[0] = mem_dims[0];
    H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

    // Select all elements in the memory buffer
    mem_start[0] = 0;
    mem_count[0] = mem_dims[0];

    H5Dread(dset_id, Hmemtype, mem_space_id, file_space_id, dxfer_plist_id, readdata);

    H5Dclose(dset_id);
    read_tend = MPI_Wtime();
    H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    free(readdata);

    // Close resources
    H5Sclose(mem_space_id);
    H5Sclose(file_space_id);
    H5Pclose(dxfer_plist_id);
    H5Tclose (Hmemtype);
    H5Fclose(file_id);
}

void create_multi_dataset(hid_t file_id, int mpi_size) {
    hid_t dcpl_ids[NUM_VARS];
    hid_t file_space_ids[NUM_VARS];
    hid_t file_ids[NUM_VARS];
    hid_t type_ids[NUM_VARS];
    hid_t dset_ids[NUM_VARS];
    hsize_t file_dims[RANK];
    hsize_t chunk_dims[RANK];

    int i;
    for (i=0; i < NUM_VARS; i++) {
        file_dims[0] = NUM_DOUBLES_PER_VAR;
        chunk_dims[0] = NUM_DOUBLES_PER_VAR / mpi_size;

        file_space_ids[i] = H5Screate_simple(1, file_dims, NULL);
        dcpl_ids[i] = H5Pcreate(H5P_DATASET_CREATE);

        H5Pset_chunk(dcpl_ids[i], 1, chunk_dims);

        H5Pset_alloc_time(dcpl_ids[i], H5D_ALLOC_TIME_MULTI);
        H5Pset_fill_time(dcpl_ids[i], H5D_FILL_TIME_NEVER);
        file_ids[i] = file_id;
        type_ids[i] = H5T_NATIVE_DOUBLE;
    }

    H5Dcreate_multi(NUM_VARS, file_ids, DATASETNAME, type_ids, file_space_ids, H5P_DEFAULT, dcpl_ids, H5P_DEFAULT, H5D_ALLOC_MULTI_ROUND_ROBIN, NULL, dset_ids);

    for (i=0; i < NUM_VARS; i++) {
        H5Sclose(file_space_ids[i]);
        H5Pclose(dcpl_ids[i]);
        H5Dclose(dset_ids[i]);
    }
}


void hacc_hdf5_with_seperate_dataset(int mpi_rank, int mpi_size, bool multi) {

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

    // 1. Create 9 datasets, one for each variable
    if (!multi) {
        dcpl_id = get_dcpl_id(mpi_size);
        file_dims[0] = BUF_SIZE_PER_VAR / sizeof(double);
        file_space_id = H5Screate_simple(1, file_dims, NULL);
        for (i = 0; i < NUM_VARS; i++) {
            dset_id = H5Dcreate2(file_id, DATASETNAME[i], H5T_NATIVE_DOUBLE, file_space_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
            H5Dclose(dset_id);
        }
        H5Sclose(file_space_id);
        H5Pclose(dcpl_id);
    } else {
        create_multi_dataset(file_id, mpi_size);
    }
    //H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    //H5Fclose(file_id);

    //file_id = open_hdf5_file(filename);

    // 2. Initialize the data for writing
    size_t NUM_DOUBLES_PER_VAR_PER_RANK = NUM_DOUBLES_PER_VAR / mpi_size;
    double *writedata = (double*) malloc(BUF_SIZE_PER_VAR * NUM_VARS / mpi_size);
    for (i = 0; i < NUM_DOUBLES_PER_VAR_PER_RANK * NUM_VARS; i++) {
        writedata[i] = (double)(mpi_rank * NUM_DOUBLES_PER_VAR_PER_RANK+ i % NUM_DOUBLES_PER_VAR_PER_RANK);
    }

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
    //print_io_stat(mpi_rank, mpi_size);
    write_tend = MPI_Wtime();
    H5Fflush(file_id, H5F_SCOPE_GLOBAL);

    MPI_Barrier(MPI_COMM_WORLD);
    free(writedata);

    // 4. read the data
    double *readdata = (double*) malloc(BUF_SIZE_PER_VAR * NUM_VARS / mpi_size);

    read_tstart = MPI_Wtime();
    for (i=0; i < NUM_VARS; i++) {
        dset_id = H5Dopen(file_id, DATASETNAME[i], H5P_DEFAULT);

        // Select column of elements in the file dataset
        file_start[0] = mpi_rank * NUM_DOUBLES_PER_VAR_PER_RANK;
        file_count[0] = NUM_DOUBLES_PER_VAR_PER_RANK;
        H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

        // Select all elements in the memory buffer
        mem_start[0] = i * NUM_DOUBLES_PER_VAR_PER_RANK;
        mem_count[0] = NUM_DOUBLES_PER_VAR_PER_RANK;
        H5Sselect_hyperslab(mem_space_id, H5S_SELECT_SET, mem_start, NULL, mem_count, NULL);

        // Read data independently
        H5Dread(dset_id, H5T_NATIVE_DOUBLE, mem_space_id, file_space_id, dxfer_plist_id, readdata);
        H5Dclose(dset_id);
    }
    read_tend = MPI_Wtime();
    H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    free(readdata);

    H5Pclose(dxfer_plist_id);
    H5Sclose(mem_space_id);
    H5Sclose(file_space_id);
    H5Fclose(file_id);
}



int main(int argc, char **argv)
{

    int mpi_size, mpi_rank;	    /* MPI variables */

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    if (argc !=3) {
        if (mpi_rank == 0)
            printf("Usange: hacc_hdf5.out [-c|-i] N(MB)\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    BUF_SIZE_PER_VAR = atoi(argv[2]) * MB;
    NUM_DOUBLES_PER_VAR = BUF_SIZE_PER_VAR / sizeof(double);

    if (strcmp(argv[1], "-i") == 0)
        hacc_hdf5_with_seperate_dataset(mpi_rank, mpi_size, false);
    else if (strcmp(argv[1], "-m") == 0)
        hacc_hdf5_with_seperate_dataset(mpi_rank, mpi_size, true);
    else
        hacc_hdf5_with_compound_type(mpi_rank, mpi_size);

    double min_tstart, max_tend, total_time, writebw;
    double max_tstart, min_tend;
    MPI_Reduce(&write_tstart, &min_tstart, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&write_tstart, &max_tstart, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&write_tend, &min_tend, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&write_tend, &max_tend, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(mpi_rank == 0) {
        printf("min tstart: %f max tstart: %f min tend: %f max tend: %f\n", min_tstart, max_tstart, min_tend, max_tend);
        total_time = max_tend - min_tstart;
        writebw = BUF_SIZE_PER_VAR*NUM_VARS/total_time/MB;
    }

    MPI_Reduce(&read_tstart, &min_tstart, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&read_tend, &max_tend, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(mpi_rank == 0) {
        total_time = max_tend - min_tstart;
        printf("Total time: %f, Aggragated Write Bandwidth: %fMB/s\tRead Bandwidth: %fMB/s\n", total_time, writebw, BUF_SIZE_PER_VAR*NUM_VARS/total_time/MB);
    }

    MPI_Finalize();
    return 0;
}
