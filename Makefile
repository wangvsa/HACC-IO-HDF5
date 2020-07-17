HDF5_DIR=$(HDF5_HOME)

CC=mpicc

all: hacc_hdf5.out hacc_mpiio.out

hacc_hdf5.out: hacc_hdf5.c
	$(CC) -o $@ $^ -I$(HDF5_DIR)/include -L$(HDF5_DIR)/lib -lhdf5

hacc_mpiio.out: hacc_mpiio.c
	$(CC) -o $@ $^

clean:
	rm -f *.o *.out
