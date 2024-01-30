#include <string.h>
#include <limits.h>
#include <mpi.h>
#include "sdecomp.h"
#include "memory.h"
#include "domain.h"
#include "array.h"
#include "fileio.h"

static int prepare(
    const domain_t * domain,
    const int nadds[NDIMS][2],
    size_t size,
    array_t * array
){
  // get total number of cells
  const size_t * mysizes = domain->mysizes;
  size_t nitems = 1;
  for(size_t dim = 0; dim < NDIMS; dim++){
    nitems *= mysizes[dim] + nadds[dim][0] + nadds[dim][1];
  }
  // for now local array size should be smaller than INT_MAX,
  //   since (4-byte) integer counters are used to sweep arrays
  //   to allow negative array indices
  int retval = 0;
  if(nitems > INT_MAX){
    printf("local array size (%zu) exceeds INT_MAX (%d)\n", nitems, INT_MAX);
    retval = 1;
  }
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_INT, MPI_MAX, comm_cart);
  if(0 != retval){
    return retval;
  }
  // assign members
  array->size = size;
  array->nadds = memory_calloc(NDIMS, 2 * sizeof(int));
  for(size_t dim = 0; dim < NDIMS; dim++){
    array->nadds[dim][0] = nadds[dim][0];
    array->nadds[dim][1] = nadds[dim][1];
  }
  array->datasize = nitems * size;
  array->data = memory_calloc(nitems, size);
  return 0;
}

static int destroy(
    array_t * array
){
  memory_free(array->nadds);
  memory_free(array->data);
  return 0;
}

// array->data, including additional (boundary & halo) cells
//   [1 - nadds[0][0] : mysizes[0] + nadds[0][1]]
//   x
//   [1 - nadds[1][0] : mysizes[1] + nadds[1][1]]
//   x
//   [1 - nadds[2][0] : mysizes[2] + nadds[2][1]]
// buf, holding only data to be written / loaded
//   [1 - nadds[0][0] : mysizes[0] + nadds[0][1]]
//   x
//   [1               : mysizes[1]              ]
//   x
//   [1               : mysizes[2]              ]

static int get_index(
    const int mysizes[NDIMS],
    const int nadds[NDIMS][2],
    const int indices[NDIMS]
){
  const int index =
    +  indices[0]
    + (indices[1] + nadds[1][0]) * (mysizes[0] + nadds[0][0] + nadds[0][1]);
  return index;
}

static int load(
    const domain_t * domain,
    const char dirname[],
    const char dsetname[],
    const char dtype[],
    array_t * array
){
  const size_t * glsizes = domain->glsizes;
  const size_t * mysizes = domain->mysizes;
  const size_t * offsets = domain->offsets;
  const int nadds[NDIMS][2] = {
    {array->nadds[0][0], array->nadds[0][1]},
    {array->nadds[1][0], array->nadds[1][1]},
  };
  const size_t size = array->size;
  char * data = array->data;
  // prepare a buffer
  char * buf = memory_calloc(
      (mysizes[0] + nadds[0][0] + nadds[0][1]) * mysizes[1],
      size
  );
  // read
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  const int retval = fileio.r_nd_parallel(
      comm_cart,
      dirname,
      dsetname,
      NDIMS,
      (int [NDIMS]){
        glsizes[1],
        glsizes[0] + nadds[0][0] + nadds[0][1],
      },
      (int [NDIMS]){
        mysizes[1],
        mysizes[0] + nadds[0][0] + nadds[0][1],
      },
      (int [NDIMS]){
        offsets[1],
        offsets[0],
      },
      dtype,
      size,
      buf
  );
  if(0 != retval){
    memory_free(buf);
    return 1;
  }
  // copy
  const int imax = mysizes[0] + nadds[0][0] + nadds[0][1];
  const int jmax = mysizes[1];
  for(int cnt = 0, j = 0; j < jmax; j++){
    for(int i = 0; i < imax; i++, cnt++){
      const int index = get_index(
          (int [NDIMS]){mysizes[0], mysizes[1]},
          nadds,
          (int [NDIMS]){i, j}
      );
      memcpy(data + size * index, buf + size * cnt, size);
    }
  }
  memory_free(buf);
  return 0;
}

static int dump(
    const domain_t * domain,
    const char dirname[],
    const char dsetname[],
    const char dtype[],
    const array_t * array
){
  const size_t * glsizes = domain->glsizes;
  const size_t * mysizes = domain->mysizes;
  const size_t * offsets = domain->offsets;
  const int nadds[NDIMS][2] = {
    {array->nadds[0][0], array->nadds[0][1]},
    {array->nadds[1][0], array->nadds[1][1]},
  };
  const size_t size = array->size;
  const char * data = array->data;
  // prepare a buffer
  char * buf = memory_calloc(
      (mysizes[0] + nadds[0][0] + nadds[0][1]) * mysizes[1],
      size
  );
  // copy
  const int imax = mysizes[0] + nadds[0][0] + nadds[0][1];
  const int jmax = mysizes[1];
  for(int cnt = 0, j = 0; j < jmax; j++){
    for(int i = 0; i < imax; i++, cnt++){
      const int index = get_index(
          (int [NDIMS]){mysizes[0], mysizes[1]},
          nadds,
          (int [NDIMS]){i, j}
      );
      memcpy(buf + size * cnt, data + size * index, size);
    }
  }
  // write
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  fileio.w_nd_parallel(
      comm_cart,
      dirname,
      dsetname,
      NDIMS,
      (int [NDIMS]){
        glsizes[1],
        glsizes[0] + nadds[0][0] + nadds[0][1],
      },
      (int [NDIMS]){
        mysizes[1],
        mysizes[0] + nadds[0][0] + nadds[0][1],
      },
      (int [NDIMS]){
        offsets[1],
        offsets[0],
      },
      dtype,
      size,
      buf
  );
  memory_free(buf);
  return 0;
}

const array_method_t array = {
  .prepare = prepare,
  .destroy = destroy,
  .load    = load,
  .dump    = dump,
};

