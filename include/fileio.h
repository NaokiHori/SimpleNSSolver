#if !defined(FILEIO_H)
#define FILEIO_H

#include <stdio.h> // FILE, size_t
#include <mpi.h>   // MPI_Datatype

typedef struct {
  // NPY datatypes, which are embedded in NPY files ("dtype" argument)
  // they are declared here and defined in src/fileio.c
  // NOTE: size (x-byte) may be wrong, depending on the architecture
  // 8-byte little-endian unsigned integer
  const char * npy_size_t;
  // 8-byte little-endian floating point
  const char * npy_double;
  // initialiser
  int (* const init)(
      void
  );
  // general-purpose file opener
  FILE * (* const fopen)(
      const char * path,
      const char * mode
  );
  // general-purpose file closer
  int (* const fclose)(
      FILE * stream
  );
  // prepare directory to be stored
  int (* const mkdir)(
      const char dirname[]
  );
  // NPY serial read (called by one process)
  int (* const r_serial)(
      const char dirname[],
      const char dsetname[],
      const size_t ndims,
      const size_t * shape,
      const char dtype[],
      const size_t size,
      void * data
  );
  // NPY serial write (called by one process)
  int (* const w_serial)(
      const char dirname[],
      const char dsetname[],
      const size_t ndims,
      const size_t * shape,
      const char dtype[],
      const size_t size,
      const void * data
  );
  // NPY parallel read of N-dimensional array (called by all processes)
  int (* const r_nd_parallel)(
      const MPI_Comm comm,
      const char dirname[],
      const char dsetname[],
      const size_t ndims,
      const int * array_of_sizes,
      const int * array_of_subsizes,
      const int * array_of_starts,
      const char dtype[],
      const size_t size,
      void * data
  );
  // NPY parallel write of N-dimensional array (called by all processes)
  int (* const w_nd_parallel)(
      const MPI_Comm comm,
      const char dirname[],
      const char dsetname[],
      const size_t ndims,
      const int * array_of_sizes,
      const int * array_of_subsizes,
      const int * array_of_starts,
      const char dtype[],
      const size_t size,
      const void * data
  );
} fileio_t;

extern const fileio_t fileio;

#endif // FILEIO_H
