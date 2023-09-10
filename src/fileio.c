#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <mpi.h>
#include "snpyio.h"
#include "memory.h"
#include "fileio.h"

// allocate and initialise string having a NPY file name
static char * create_npyfname(
    const char dirname[],
    const char dsetname[]
){
  if(NULL == dirname){
    fprintf(stderr, "dirname is NULL\n");
    return NULL;
  }
  if(NULL == dsetname){
    fprintf(stderr, "dsetname is NULL\n");
    return NULL;
  }
  const char slash[] = {"/"};
  const char suffix[] = {".npy"};
  const size_t nchars =
    + strlen( dirname)
    + strlen(   slash)
    + strlen(dsetname)
    + strlen(  suffix);
  char * fname = memory_calloc(nchars + 2, sizeof(char));
  snprintf(fname, nchars + 1, "%s%s%s%s", dirname, slash, dsetname, suffix);
  fname[nchars + 1] = '\0';
  return fname;
}

static int mpi_file_open(
    const MPI_Comm comm,
    char * fname,
    int amode,
    MPI_File * fh
){
  const int mpi_error_code = MPI_File_open(comm, fname, amode, MPI_INFO_NULL, fh);
  if(MPI_SUCCESS != mpi_error_code){
    char string[MPI_MAX_ERROR_STRING] = {'\0'};
    int resultlen = 0;
    MPI_Error_string(mpi_error_code, string, &resultlen);
    fprintf(stderr, "%s: %s\n", fname, string);
    memory_free(fname);
    return 1;
  }
  return 0;
}

static int get_count(
    const size_t ndims,
    const int * mysizes
){
  int count = 1;
  for(size_t n = 0; n < ndims; n++){
    count *= mysizes[n];
  }
  return count;
}

static int prepare_view(
    const size_t ndims,
    const int * glsizes,
    const int * mysizes,
    const int * offsets,
    MPI_File fh,
    const size_t header_size,
    const MPI_Datatype basetype,
    MPI_Datatype * filetype
){
  // create data type and set file view
  MPI_Type_create_subarray((int)ndims, glsizes, mysizes, offsets, MPI_ORDER_C, basetype, filetype);
  MPI_Type_commit(filetype);
  MPI_File_set_view(fh, (MPI_Offset)header_size, basetype, *filetype, "native", MPI_INFO_NULL);
  return 0;
}

static int destroy_view(
    MPI_Datatype * filetype
){
  // clean-up datatype
  MPI_Type_free(filetype);
  return 0;
}

static int init(
    void
){
  const int root = 0;
  int myrank = root;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  const size_t sizeof_size_t = sizeof(size_t);
  const size_t sizeof_double = sizeof(double);
  FILE * stream = stderr;
  if(8 != sizeof_size_t){
    if(root == myrank){
      fprintf(stream, "NPY data type '<u8' and sizeof(size_t): %zu mismatch\n", sizeof_size_t);
      fflush(stream);
    }
    return 1;
  }
  if(8 != sizeof_double){
    if(root == myrank){
      fprintf(stream, "NPY data type '<i4' and sizeof(double): %zu mismatch\n", sizeof_double);
      fflush(stream);
    }
    return 1;
  }
  return 0;
}

static FILE * fopen_(
    const char * path,
    const char * mode
){
  errno = 0;
  FILE * stream = fopen(path, mode);
  if(NULL == stream){
    perror(path);
  }
  return stream;
}

static int fclose_(
    FILE * stream
){
  if(NULL == stream){
    return 1;
  }
  fclose(stream);
  return 0;
}

// create directory
static int mkdir_(
    const char dirname[]
){
  // NOTE: call this function ONLY from the main process
  // NOTE: continue even if failed,
  //   since we want to override previous data (errorcode: EEXIST)
  // RWX masks for user, group, and others (0o777, ref. "man 2 chmod")
  const mode_t mode = S_IRWXU | S_IRWXG | S_IRWXO;
  if(0 != mkdir(dirname, mode)){
    perror(dirname);
    return 1;
  }
  return 0;
}

// wrapper function of snpyio_r_header with error handling
static size_t r_npy_header(
    const char fname[],
    const size_t ndims,
    const size_t * shape,
    const char * dtype,
    const bool is_fortran_order
){
  const char msg[] = {"NPY header read failed"};
  FILE * fp = fopen_(fname, "r");
  if(NULL == fp){
    return 0;
  }
  // load header, return header size when succeeded, return 0 otherwise
  size_t ndims_ = 0;
  size_t * shape_ = NULL;
  char * dtype_ = NULL;
  bool is_fortran_order_ = false;
  size_t header_size = snpyio_r_header(&ndims_, &shape_, &dtype_, &is_fortran_order_, fp);
  fclose_(fp);
  // check arguments, return loaded header size when all OK, return 0 otherwise
  // ndims
  if(ndims != ndims_){
    fprintf(stderr, "%s(%s), ndims: %zu expected, %zu obtained\n", msg, fname, ndims, ndims_);
    header_size = 0;
    goto err_hndl;
  }
  // shape (for each dimension)
  for(size_t n = 0; n < ndims; n++){
    if(shape[n] != shape_[n]){
      fprintf(stderr, "%s(%s), shape[%zu]: %zu expected, %zu obtained\n", msg, fname, n, shape[n], shape_[n]);
      header_size = 0;
      goto err_hndl;
    }
  }
  // dtype
  if(0 != strcmp(dtype, dtype_)){
    fprintf(stderr, "%s(%s), dtype: %s expected, %s obtained\n", msg, fname, dtype, dtype_);
    header_size = 0;
    goto err_hndl;
  }
  // is_fortran_order
  if(is_fortran_order != is_fortran_order_){
    fprintf(stderr, "%s(%s), is_fortran_order: %s expected, %s obtained\n", msg, fname, is_fortran_order ? "true" : "false", is_fortran_order_ ? "true" : "false");
    header_size = 0;
    goto err_hndl;
  }
err_hndl:
  memory_free(shape_);
  memory_free(dtype_);
  return header_size;
}

// wrapper function of snpyio_w_header with error handling
static size_t w_npy_header(
    const char fname[],
    const size_t ndims,
    const size_t * shape,
    const char dtype[],
    const bool is_fortran_order
){
  const char msg[] = {"NPY header write failed"};
  FILE * fp = fopen_(fname, "w");
  if(NULL == fp){
    return 0;
  }
  const size_t header_size = snpyio_w_header(ndims, shape, dtype, is_fortran_order, fp);
  if(0 == header_size){
    fprintf(stderr, "%s(%s)\n", msg, fname);
  }
  fclose_(fp);
  return header_size;
}

/**
 * @brief read data from a npy file, by one process
 * @param[in]  dirname  : name of directory in which a target npy file is contained
 * @param[in]  dsetname : name of dataset
 * @param[in]  ndims    : number of dimensions of dataset
 * @param[in]  shape    : shape of dataset
 * @param[in]  dtype    : datatype, e.g. '<f8'
 * @param[in]  size     : size of each element
 * @param[out] data     : pointer to the data to be loaded
 */
static int r_serial(
    const char dirname[],
    const char dsetname[],
    const size_t ndims,
    const size_t * shape,
    const char dtype[],
    const size_t size,
    void * data
){
  char * fname = create_npyfname(dirname, dsetname);
  const size_t header_size = r_npy_header(fname, ndims, shape, dtype, false);
  if(0 == header_size){
    fprintf(stderr, "%s: NPY header load failed\n", fname);
    memory_free(fname);
    return 1;
  }
  FILE * fp = fopen_(fname, "r");
  if(NULL == fp){
    memory_free(fname);
    return 1;
  }
  if(0 != fseek(fp, (long)header_size, SEEK_SET)){
    fprintf(stderr, "%s: fseek failed\n", fname);
    fclose_(fp);
    memory_free(fname);
    return 1;
  }
  size_t nitems = 1;
  for(size_t dim = 0; dim < ndims; dim++){
    nitems *= shape[dim];
  }
  const size_t nitems_ = fread(data, size, nitems, fp);
  if(nitems_ != nitems){
    fprintf(stderr, "%s: fread failed (%zu vs %zu)\n", fname, nitems_, nitems);
    fclose_(fp);
    memory_free(fname);
    return 1;
  }
  fclose_(fp);
  memory_free(fname);
  return 0;
}

/**
 * @brief write data to a npy file, by one process
 * @param[in] dirname  : name of directory in which a target npy file is contained
 * @param[in] dsetname : name of dataset
 * @param[in] ndims    : number of dimensions of dataset
 * @param[in] shape    : shape of dataset
 * @param[in] dtype    : datatype, e.g. '<f8'
 * @param[in] size     : size of each element
 * @param[in] data     : pointer to the data to be written
 */
static int w_serial(
    const char dirname[],
    const char dsetname[],
    const size_t ndims,
    const size_t * shape,
    const char dtype[],
    const size_t size,
    const void * data
){
  char * fname = create_npyfname(dirname, dsetname);
  const size_t header_size = w_npy_header(fname, ndims, shape, dtype, false);
  if(0 == header_size){
    memory_free(fname);
    return 1;
  }
  FILE * fp = fopen_(fname, "a");
  if(NULL == fp){
    memory_free(fname);
    return 1;
  }
  if(0 != fseek(fp, (long)header_size, SEEK_SET)){
    fprintf(stderr, "%s: fseek failed\n", fname);
    fclose_(fp);
    memory_free(fname);
    return 1;
  }
  size_t nitems = 1;
  for(size_t dim = 0; dim < ndims; dim++){
    nitems *= shape[dim];
  }
  const size_t nitems_ = fwrite(data, size, nitems, fp);
  if(nitems_ != nitems){
    fprintf(stderr, "%s: fwrite failed\n", fname);
    fclose_(fp);
    memory_free(fname);
    return 1;
  }
  fclose_(fp);
  memory_free(fname);
  return 0;
}

/**
 * @brief read N-dimensional data from a npy file, by all processes
 * @param[in]  comm     : communicator to which all processes calling this function belong
 * @param[in]  dirname  : name of directory in which a target npy file is contained
 * @param[in]  dsetname : name of dataset
 * @param[in]  ndims    : number of dimensions of the array
 * @param[in]  glsizes  : global sizes   of the dataset
 * @param[in]  mysizes  : local  sizes   of the dataset
 * @param[in]  offsets  : local  offsets of the dataset
 * @param[in]  dtype    : NPY data type
 * @param[in]  size     : size of each element
 * @param[out] data     : pointer to the data to be loaded
 */
static int r_nd_parallel(
    const MPI_Comm comm,
    const char dirname[],
    const char dsetname[],
    const size_t ndims,
    const int * glsizes,
    const int * mysizes,
    const int * offsets,
    const char dtype[],
    const size_t size,
    void * data
){
  const int root = 0;
  int myrank = root;
  MPI_Comm_rank(comm, &myrank);
  char * fname = create_npyfname(dirname, dsetname);
  // check header by main process
  size_t header_size = 0;
  if(root == myrank){
    // set values which are expected to be in NPY file
    size_t * shape = memory_calloc(ndims, sizeof(size_t));
    for(size_t dim = 0; dim < ndims; dim++){
      shape[dim] = (size_t)glsizes[dim];
    }
    header_size = r_npy_header(fname, ndims, shape, dtype, false);
    memory_free(shape);
  }
  // share result
  MPI_Bcast(&header_size, sizeof(size_t), MPI_BYTE, root, comm);
  if(0 == header_size){
    memory_free(fname);
    return 1;
  }
  // open file
  MPI_File fh = NULL;
  if(0 != mpi_file_open(comm, fname, MPI_MODE_RDONLY, &fh)) return 1;
  // prepare file view
  MPI_Datatype basetype = MPI_BYTE;
  MPI_Datatype filetype = MPI_DATATYPE_NULL;
  MPI_Type_contiguous(size, basetype, &basetype);
  MPI_Type_commit(&basetype);
  prepare_view((int)ndims, glsizes, mysizes, offsets, fh, header_size, basetype, &filetype);
  // get number of elements which are locally read
  const int count = get_count(ndims, mysizes);
  // read
  MPI_File_read_all(fh, data, count, basetype, MPI_STATUS_IGNORE);
  // clean-up file view
  MPI_Type_free(&basetype);
  destroy_view(&filetype);
  // close file
  MPI_File_close(&fh);
  memory_free(fname);
  return 0;
}

/**
 * @brief write N-dimensional data to a npy file, by all processes
 * @param[in] comm     : communicator to which all processes calling this function belong
 * @param[in] dirname  : name of directory in which a target npy file is contained
 * @param[in] dsetname : name of dataset
 * @param[in] ndims    : number of dimensions of the array
 * @param[in] glsizes  : global sizes   of the dataset
 * @param[in] mysizes  : local  sizes   of the dataset
 * @param[in] offsets  : local  offsets of the dataset
 * @param[in] dtype    : NPY data type
 * @param[in] size     : size of each element
 * @param[in] data     : pointer to the data to be written
 */
static int w_nd_parallel(
    const MPI_Comm comm,
    const char dirname[],
    const char dsetname[],
    const size_t ndims,
    const int * glsizes,
    const int * mysizes,
    const int * offsets,
    const char dtype[],
    const size_t size,
    const void * data
){
  const int root = 0;
  int myrank = root;
  MPI_Comm_rank(comm, &myrank);
  char * fname = create_npyfname(dirname, dsetname);
  // check header by main process
  size_t header_size = 0;
  if(root == myrank){
    // set values which are expected to be in NPY file
    size_t * shape = memory_calloc(ndims, sizeof(size_t));
    for(size_t dim = 0; dim < ndims; dim++){
      shape[dim] = (size_t)glsizes[dim];
    }
    header_size = w_npy_header(fname, ndims, shape, dtype, false);
    memory_free(shape);
  }
  // share result
  MPI_Bcast(&header_size, sizeof(size_t), MPI_BYTE, root, comm);
  if(0 == header_size){
    memory_free(fname);
    return 1;
  }
  // open file
  MPI_File fh = NULL;
  if(0 != mpi_file_open(comm, fname, MPI_MODE_CREATE | MPI_MODE_RDWR, &fh)) return 1;
  // prepare file view
  MPI_Datatype basetype = MPI_BYTE;
  MPI_Datatype filetype = MPI_DATATYPE_NULL;
  MPI_Type_contiguous(size, basetype, &basetype);
  MPI_Type_commit(&basetype);
  prepare_view((int)ndims, glsizes, mysizes, offsets, fh, header_size, basetype, &filetype);
  // get number of elements which are locally written
  const int count = get_count(ndims, mysizes);
  // write
  MPI_File_write_all(fh, data, count, basetype, MPI_STATUS_IGNORE);
  // clean-up file view
  MPI_Type_free(&basetype);
  destroy_view(&filetype);
  // close file
  MPI_File_close(&fh);
  memory_free(fname);
  return 0;
}

const fileio_t fileio = {
  .npy_size_t = "'<u8'",
  .npy_double = "'<f8'",
  .init = init,
  .fopen = fopen_,
  .fclose = fclose_,
  .mkdir = mkdir_,
  .r_serial = r_serial,
  .w_serial = w_serial,
  .r_nd_parallel = r_nd_parallel,
  .w_nd_parallel = w_nd_parallel,
};

