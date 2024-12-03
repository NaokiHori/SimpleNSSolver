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

#define REPORT_ERROR(...) \
  do { \
    fprintf(stderr, "[%s:%d] ", __FILE__, __LINE__); \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, "\n"); \
  } while (0)

static const char NPY_SIZE_T[] = "'<u8'";
static const char NPY_DOUBLE[] = "'<f8'";

static char * create_npy_file_name(
    const char directory_name[],
    const char dataset_name[]
) {
  if (NULL == directory_name) {
    REPORT_ERROR("directory_name is NULL");
    return NULL;
  }
  if (NULL == dataset_name) {
    REPORT_ERROR("dataset_name is NULL");
    return NULL;
  }
  // avoid a double-slash
  const char * const slash = (directory_name[strlen(directory_name) - 1] == '/') ? "" : "/";
  const char suffix[] = ".npy";
  const size_t nchars =
    + strlen(directory_name)
    + strlen(slash)
    + strlen(dataset_name)
    + strlen(suffix);
  // a nul character should be added to the tail: +1
  char * const file_name = memory_calloc(nchars + 1, sizeof(char));
  const int error_code = snprintf(file_name, nchars + 1, "%s%s%s%s", directory_name, slash, dataset_name, suffix);
  if (error_code < 0 || nchars + 1 <= (size_t)error_code) {
      REPORT_ERROR("snprintf failed or output truncated");
      free(file_name);
      return NULL;
  }
  return file_name;
}

static int mpi_file_open(
    const MPI_Comm comm,
    char * const file_name,
    int amode,
    MPI_File * const fh
) {
  const int error_code = MPI_File_open(comm, file_name, amode, MPI_INFO_NULL, fh);
  if (MPI_SUCCESS != error_code) {
    char error_string[MPI_MAX_ERROR_STRING + 1] = {'\0'};
    int string_length = 0;
    MPI_Error_string(error_code, error_string, &string_length);
    error_string[string_length] = '\0';
    REPORT_ERROR("%s", error_string);
    return 1;
  }
  return 0;
}

static int get_count(
    const size_t ndims,
    const int * const mysizes
) {
  int count = 1;
  for (size_t n = 0; n < ndims; n++) {
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
) {
  // create data type and set file view
  MPI_Type_create_subarray((int)ndims, glsizes, mysizes, offsets, MPI_ORDER_C, basetype, filetype);
  MPI_Type_commit(filetype);
  MPI_File_set_view(fh, (MPI_Offset)header_size, basetype, *filetype, "native", MPI_INFO_NULL);
  return 0;
}

static int destroy_view(
    MPI_Datatype * const filetype
) {
  // clean-up datatype
  MPI_Type_free(filetype);
  return 0;
}

static int init(
    void
) {
  const int root = 0;
  int myrank = root;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  const size_t sizeof_size_t = sizeof(size_t);
  const size_t sizeof_double = sizeof(double);
  if (8 != sizeof_size_t) {
    if (root == myrank) {
      REPORT_ERROR("NPY data type %s and sizeof(size_t): %zu mismatch\n", NPY_SIZE_T, sizeof_size_t);
    }
    return 1;
  }
  if (8 != sizeof_double) {
    if (root == myrank) {
      REPORT_ERROR("NPY data type %s and sizeof(double): %zu mismatch\n", NPY_DOUBLE, sizeof_double);
    }
    return 1;
  }
  return 0;
}

static FILE * fopen_(
    const char * const path,
    const char * const mode
) {
  errno = 0;
  FILE * const stream = fopen(path, mode);
  if (NULL == stream) {
    perror(path);
  }
  return stream;
}

static int fclose_(
    FILE * const stream
) {
  if (NULL == stream) {
    return 1;
  }
  fclose(stream);
  return 0;
}

// create directory
static int mkdir_(
    const char directory_name[]
) {
  // NOTE: call this function ONLY from the main process
  // NOTE: continue even if failed,
  //         since we want to override previous data (errorcode: EEXIST)
  // RWX masks for user, group, and others (0o777, ref. "man 2 chmod")
  const mode_t mode = S_IRWXU | S_IRWXG | S_IRWXO;
  errno = 0;
  if (0 != mkdir(directory_name, mode)) {
    perror(directory_name);
    return 1;
  }
  return 0;
}

// wrapper function of snpyio_r_header with error handling
static int r_npy_header(
    const char fname[],
    const size_t ndims,
    const size_t * shape,
    const char * dtype,
    const bool is_fortran_order,
    size_t * header_size
) {
  int error_code = 0;
  const char msg[] = {"NPY header read failed"};
  errno = 0;
  FILE * fp = fopen_(fname, "r");
  if (NULL == fp) {
    perror(fname);
    return 1;
  }
  // load header, return header size when succeeded, return 0 otherwise
  size_t ndims_ = 0;
  size_t * shape_ = NULL;
  char * dtype_ = NULL;
  bool is_fortran_order_ = false;
  *header_size = 0;
  error_code = snpyio_r_header(&ndims_, &shape_, &dtype_, &is_fortran_order_, fp, header_size);
  fclose_(fp);
  if (0 != error_code) {
    REPORT_ERROR("%s(%s), snpyio_r_header failed\n", msg, fname);
    goto err_hndl;
  }
  // check arguments, return loaded header size when all OK, return 0 otherwise
  // ndims
  if (ndims != ndims_) {
    REPORT_ERROR("%s(%s), ndims: %zu expected, %zu obtained\n", msg, fname, ndims, ndims_);
    error_code = 1;
    goto err_hndl;
  }
  // shape (for each dimension)
  for(size_t n = 0; n < ndims; n++) {
    if (shape[n] != shape_[n]) {
      REPORT_ERROR("%s(%s), shape[%zu]: %zu expected, %zu obtained\n", msg, fname, n, shape[n], shape_[n]);
      error_code = 1;
      goto err_hndl;
    }
  }
  // dtype
  if (0 != strcmp(dtype, dtype_)) {
    REPORT_ERROR("%s(%s), dtype: %s expected, %s obtained\n", msg, fname, dtype, dtype_);
    error_code = 1;
    goto err_hndl;
  }
  // is_fortran_order
  if (is_fortran_order != is_fortran_order_) {
    REPORT_ERROR("%s(%s), is_fortran_order: %s expected, %s obtained\n", msg, fname, is_fortran_order ? "true" : "false", is_fortran_order_ ? "true" : "false");
    error_code = 1;
    goto err_hndl;
  }
err_hndl:
  memory_free(shape_);
  memory_free(dtype_);
  return error_code;
}

// wrapper function of snpyio_w_header with error handling
static int w_npy_header(
    const char fname[],
    const size_t ndims,
    const size_t * shape,
    const char dtype[],
    const bool is_fortran_order,
    size_t * header_size
) {
  int error_code = 0;
  const char msg[] = {"NPY header write failed"};
  errno = 0;
  FILE * fp = fopen_(fname, "w");
  if (NULL == fp) {
    perror(fname);
    return 1;
  }
  *header_size = 0;
  error_code = snpyio_w_header(ndims, shape, dtype, is_fortran_order, fp, header_size);
  if (0 != error_code) {
    REPORT_ERROR("%s(%s) snpyio_w_header failed\n", msg, fname);
  }
  fclose_(fp);
  return error_code;
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
) {
  char * fname = create_npy_file_name(dirname, dsetname);
  size_t header_size = 0;
  if (0 != r_npy_header(fname, ndims, shape, dtype, false, &header_size)) {
    REPORT_ERROR("%s: NPY header load failed\n", fname);
    memory_free(fname);
    return 1;
  }
  FILE * fp = fopen_(fname, "r");
  if (NULL == fp) {
    memory_free(fname);
    return 1;
  }
  if (0 != fseek(fp, (long)header_size, SEEK_SET)) {
    REPORT_ERROR("%s: fseek failed\n", fname);
    fclose_(fp);
    memory_free(fname);
    return 1;
  }
  size_t nitems = 1;
  for(size_t dim = 0; dim < ndims; dim++) {
    nitems *= shape[dim];
  }
  const size_t nitems_ = fread(data, size, nitems, fp);
  if (nitems_ != nitems) {
    REPORT_ERROR("%s: fread failed (%zu vs %zu)\n", fname, nitems_, nitems);
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
) {
  char * fname = create_npy_file_name(dirname, dsetname);
  size_t header_size = 0;
  if (0 != w_npy_header(fname, ndims, shape, dtype, false, &header_size)) {
    memory_free(fname);
    return 1;
  }
  FILE * fp = fopen_(fname, "a");
  if (NULL == fp) {
    memory_free(fname);
    return 1;
  }
  if (0 != fseek(fp, (long)header_size, SEEK_SET)) {
    REPORT_ERROR("%s: fseek failed\n", fname);
    fclose_(fp);
    memory_free(fname);
    return 1;
  }
  size_t nitems = 1;
  for(size_t dim = 0; dim < ndims; dim++) {
    nitems *= shape[dim];
  }
  const size_t nitems_ = fwrite(data, size, nitems, fp);
  if (nitems_ != nitems) {
    REPORT_ERROR("%s: fwrite failed\n", fname);
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
) {
  int error_code = 0;
  const int root = 0;
  int myrank = root;
  MPI_Comm_rank(comm, &myrank);
  char * fname = create_npy_file_name(dirname, dsetname);
  // check header by main process
  size_t header_size = 0;
  if (root == myrank) {
    // set values which are expected to be in NPY file
    size_t * shape = memory_calloc(ndims, sizeof(size_t));
    for(size_t dim = 0; dim < ndims; dim++) {
      shape[dim] = (size_t)glsizes[dim];
    }
    error_code = r_npy_header(fname, ndims, shape, dtype, false, &header_size);
    memory_free(shape);
  }
  // share result
  MPI_Bcast(&error_code, sizeof(int), MPI_BYTE, root, comm);
  MPI_Bcast(&header_size, sizeof(size_t), MPI_BYTE, root, comm);
  if (0 != error_code) {
    goto err_hndl;
  }
  // open file
  MPI_File fh = NULL;
  if (0 != mpi_file_open(comm, fname, MPI_MODE_RDONLY, &fh)) {
    error_code = 1;
    REPORT_ERROR("mpi_file_open failed");
    goto err_hndl;
  }
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
err_hndl:
  memory_free(fname);
  return error_code;
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
) {
  int error_code = 0;
  const int root = 0;
  int myrank = root;
  MPI_Comm_rank(comm, &myrank);
  char * fname = create_npy_file_name(dirname, dsetname);
  // check header by main process
  size_t header_size = 0;
  if (root == myrank) {
    // set values which are expected to be in NPY file
    size_t * shape = memory_calloc(ndims, sizeof(size_t));
    for(size_t dim = 0; dim < ndims; dim++) {
      shape[dim] = (size_t)glsizes[dim];
    }
    error_code = w_npy_header(fname, ndims, shape, dtype, false, &header_size);
    memory_free(shape);
  }
  // share result
  MPI_Bcast(&error_code, sizeof(int), MPI_BYTE, root, comm);
  MPI_Bcast(&header_size, sizeof(size_t), MPI_BYTE, root, comm);
  if (0 != error_code) {
    goto err_hndl;
  }
  // open file
  MPI_File fh = NULL;
  if (0 != mpi_file_open(comm, fname, MPI_MODE_CREATE | MPI_MODE_RDWR, &fh)) {
    goto err_hndl;
    return 1;
  }
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
err_hndl:
  memory_free(fname);
  return error_code;
}

const fileio_t fileio = {
  .npy_size_t = NPY_SIZE_T,
  .npy_double = NPY_DOUBLE,
  .init = init,
  .fopen = fopen_,
  .fclose = fclose_,
  .mkdir = mkdir_,
  .r_serial = r_serial,
  .w_serial = w_serial,
  .r_nd_parallel = r_nd_parallel,
  .w_nd_parallel = w_nd_parallel,
};

