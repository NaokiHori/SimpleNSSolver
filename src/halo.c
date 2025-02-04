#include <stdio.h>
#include <mpi.h>
#include "array.h"
#include "domain.h"
#include "halo.h"

// fixed parameters
// since data type is defined, number of items is 1
static const int nitems = 1;
// same tag can be used since I use blocking communication
static const int tag = 0;

// assume the given data type has not been initialised yet
static const MPI_Datatype dtype_default = MPI_DOUBLE;

// communicate halo cells with the y-neighbour processes
// NOTE: send boundary cells for simplicity
int halo_communicate_in_y(
    const domain_t * domain,
    MPI_Datatype * dtype,
    array_t * array
){
  // extract communicator
  const sdecomp_info_t * info = domain->info;
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(info, &comm_cart);
  // check negative / positive neighbour ranks
  int neighbours[2] = {MPI_PROC_NULL, MPI_PROC_NULL};
  sdecomp.get_neighbours(info, SDECOMP_X1PENCIL, SDECOMP_YDIR, neighbours);
  // array size (with halo and boundary cells)
  const int isize_ = domain->mysizes[0] + array->nadds[0][0] + array->nadds[0][1];
  const int jsize_ = domain->mysizes[1] + array->nadds[1][0] + array->nadds[1][1];
  // number of halo cells
  // this function assumes same number of halo cells
  //   in the negative / positive directions
  if(array->nadds[1][0] != array->nadds[1][1]){
    printf("%s: number of halo cells in y (%d and %d) mismatch\n",
        __func__, array->nadds[1][0], array->nadds[1][1]);
    return 1;
  }
  const int nhalos_y = array->nadds[1][0];
  // define datatype in y
  if(dtype_default == *dtype){
    MPI_Type_contiguous(
        isize_ * nhalos_y,
        *dtype,
        dtype
    );
    MPI_Type_commit(dtype);
  }
  // send to positive, receive from negative
  {
    const int sindices[NDIMS] = {0, jsize_ - 2 * nhalos_y};
    const int rindices[NDIMS] = {0,          0 * nhalos_y};
    const size_t soffset = sindices[0] + isize_ * sindices[1];
    const size_t roffset = rindices[0] + isize_ * rindices[1];
    MPI_Sendrecv(
      (char *)array->data + array->size * soffset, nitems, *dtype, neighbours[1], tag,
      (char *)array->data + array->size * roffset, nitems, *dtype, neighbours[0], tag,
      comm_cart, MPI_STATUS_IGNORE
    );
  }
  // send to negative, receive from positive
  {
    const int sindices[NDIMS] = {0,          1 * nhalos_y};
    const int rindices[NDIMS] = {0, jsize_ - 1 * nhalos_y};
    const size_t soffset = sindices[0] + isize_ * sindices[1];
    const size_t roffset = rindices[0] + isize_ * rindices[1];
    MPI_Sendrecv(
      (char *)array->data + array->size * soffset, nitems, *dtype, neighbours[0], tag,
      (char *)array->data + array->size * roffset, nitems, *dtype, neighbours[1], tag,
      comm_cart, MPI_STATUS_IGNORE
    );
  }
  return 0;
}

