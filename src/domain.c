#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include "sdecomp.h"
#include "memory.h"
#include "domain.h"
#include "fileio.h"
#include "array_macros/domain/xf.h"
#include "array_macros/domain/xc.h"
#include "array_macros/domain/dxf.h"
#include "array_macros/domain/dxc.h"

/**
 * @brief load members in domain_t
 * @param[in]  dirname : name of directory from which data is loaded
 * @param[out] domain  : global domain sizes and resolutions
 * @return             : error code
 */
static int domain_load(
    const char dirname[],
    domain_t * domain
){
  size_t * glsizes = domain->glsizes;
  double * restrict lengths = domain->lengths;
  double * restrict * restrict xf = &domain->xf;
  double * restrict * restrict xc = &domain->xc;
  if(0 != fileio.r_serial(dirname, "glsizes", 1, (size_t [1]){NDIMS}, fileio.npy_size_t, sizeof(size_t), glsizes)){
    return 1;
  }
  if(0 != fileio.r_serial(dirname, "lengths", 1, (size_t [1]){NDIMS}, fileio.npy_double, sizeof(double), lengths)){
    return 1;
  }
  *xf = memory_calloc(glsizes[0] + 1, sizeof(double));
  if(0 != fileio.r_serial(dirname, "xf", 1, (size_t [1]){glsizes[0] + 1}, fileio.npy_double, sizeof(double), *xf)){
    return 1;
  }
  *xc = memory_calloc(glsizes[0] + 2, sizeof(double));
  if(0 != fileio.r_serial(dirname, "xc", 1, (size_t [1]){glsizes[0] + 2}, fileio.npy_double, sizeof(double), *xc)){
    return 1;
  }
  return 0;
}

/**
 * @brief save members in domain_t
 * @param[in] dirname : name of directory to which data is saved
 * @param[in] domain  : global domain sizes and resolutions
 * @return            : error code
 */
int domain_save(
    const char dirname[],
    const domain_t * domain
){
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  // since this is a serial operation,
  //   other processes are not involved
  if(root != myrank){
    return 0;
  }
  const size_t * glsizes = domain->glsizes;
  const double * lengths = domain->lengths;
  const double * xf      = domain->xf;
  const double * xc      = domain->xc;
  fileio.w_serial(dirname, "glsizes", 1, (size_t [1]){NDIMS}, fileio.npy_size_t, sizeof(size_t), glsizes);
  fileio.w_serial(dirname, "lengths", 1, (size_t [1]){NDIMS}, fileio.npy_double, sizeof(double), lengths);
  fileio.w_serial(dirname, "xf", 1, (size_t [1]){glsizes[0] + 1}, fileio.npy_double, sizeof(double), xf);
  fileio.w_serial(dirname, "xc", 1, (size_t [1]){glsizes[0] + 2}, fileio.npy_double, sizeof(double), xc);
  return 0;
}

int domain_check_x_grid_is_uniform(
    const domain_t * domain,
    bool * x_grid_is_uniform
){
  static bool is_checked = false;
  static bool is_uniform = false;
  if(!is_checked){
    const size_t isize = domain->glsizes[0];
    const double * dxf = domain->dxf;
    double extrema[2] = {+1. * DBL_MAX, -1. * DBL_MAX};
    for(size_t i = 0; i < isize; i++){
      extrema[0] = fmin(extrema[0], dxf[i]);
      extrema[1] = fmax(extrema[1], dxf[i]);
    }
    if(fabs(extrema[1] - extrema[0]) < 1.e-15){
      is_uniform = true;
    }else{
      is_uniform = false;
    }
    is_checked = true;
  }
  *x_grid_is_uniform = is_uniform;
  return 0;
}

/**
 * @brief define face-to-face distances in x direction
 * @param[in] isize : number of cell-centers in x direction (boundary excluded)
 * @param[in] xf    : cell-face positions in x direction
 * @return          : face-to-face distances in x direction
 */
static double * allocate_and_init_dxf(
    const int isize,
    const double * xf
){
  // dxf: distance from cell face to cell face
  // NOTE: since xf has "isize + 1" items,
  //   dxf, which tells the distance of the two neighbouring cell faces,
  //   has "isize" elements, whose index starts from 1
  const size_t nitems = isize;
  double * dxf = memory_calloc(nitems, sizeof(double));
  for(size_t i = 1; i <= nitems; i++){
    DXF(i  ) = XF(i+1) - XF(i  );
  }
  return dxf;
}

/**
 * @brief define center-to-center distances in x direction
 * @param[in] isize : number of cell-centers in x direction (boundary excluded)
 * @param[in] xc    : cell-center positions in x direction
 * @return          : center-to-center distances in x direction
 */
static double * allocate_and_init_dxc(
    const int isize,
    const double * xc
){
  // dxc: distance from cell center to cell center (generally)
  // NOTE: since xc has "isize + 2" items,
  //   dxc, which tells the distance of the two neighbouring cell centers,
  //   has "isize + 1" elements, whose index starts from 1
  const size_t nitems = isize + 1;
  double * dxc = memory_calloc(nitems, sizeof(double));
  for(size_t i = 1; i <= nitems; i++){
    DXC(i  ) = XC(i  ) - XC(i-1);
  }
  return dxc;
}

static void report(
    const domain_t * domain
){
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(root == myrank){
    printf("DOMAIN\n");
    for(sdecomp_dir_t dim = 0; dim < NDIMS; dim++){
      printf("\tglsizes[%u]: %zu\n", dim, domain->glsizes[dim]);
    }
    for(sdecomp_dir_t dim = 0; dim < NDIMS; dim++){
      printf("\tlengths[%u]: % .7e\n", dim, domain->lengths[dim]);
    }
    fflush(stdout);
  }
}

/**
 * @brief constructor of the structure
 * @param[in]  dirname_ic : name of directory in which initial conditions are stored
 * @param[out] domain     : structure being allocated and initalised
 * @return                : (success) 0
 *                          (failure) non-zero value
 */
int domain_init(
    const char dirname_ic[],
    domain_t * domain
){
  sdecomp_info_t ** info    = &domain->info;
  size_t * restrict glsizes =  domain->glsizes;
  size_t * restrict mysizes =  domain->mysizes;
  size_t * restrict offsets =  domain->offsets;
  double * restrict lengths =  domain->lengths;
  double * restrict * xf    = &domain->xf;
  double * restrict * xc    = &domain->xc;
  double * restrict * dxf   = &domain->dxf;
  double * restrict * dxc   = &domain->dxc;
  double * restrict   dy    = &domain->dy;
  // load spatial information
  if(0 != domain_load(dirname_ic, domain)){
    return 1;
  }
  // compute grid sizes
  // allocate and initialise x coordinates
  *dxf = allocate_and_init_dxf(glsizes[0], *xf);
  *dxc = allocate_and_init_dxc(glsizes[0], *xc);
  // grid sizes in homogeneous directions
  *dy = lengths[1] / glsizes[1];
  // initialise sdecomp to distribute the domain
  if(0 != sdecomp.construct(
        MPI_COMM_WORLD,
        NDIMS,
        (size_t [NDIMS]){0, 0},
        (bool [NDIMS]){false, true},
        info
  )) return 1;
  // local array sizes and offsets
  for(size_t dim = 0; dim < NDIMS; dim++){
    sdecomp.get_pencil_mysize(*info, SDECOMP_X1PENCIL, dim, glsizes[dim], mysizes + dim);
    sdecomp.get_pencil_offset(*info, SDECOMP_X1PENCIL, dim, glsizes[dim], offsets + dim);
  }
  report(domain);
  return 0;
}

