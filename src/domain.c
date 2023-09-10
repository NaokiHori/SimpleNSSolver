#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#if NDIMS == 3
#include <fftw3.h>
#include "timer.h"
#endif
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

#if NDIMS == 3
static int get_ndigits(
    int num
){
  // just for pretty print
  // e.g. num =    3 -> return 1
  // e.g. num =   13 -> return 2
  // e.g. num = 1234 -> return 4
  if(num < 0){
    return 0;
  }
  int retval = 1;
  while(num /= 10){
    retval++;
  }
  return retval;
}

// optimise MPI domain decomposition,
//   i.e. minimise all-to-all time
static int optimise_sdecomp_init(
    const size_t * glsizes,
    sdecomp_info_t ** info_optimum
){
  // periodicity in each dimension,
  //   which are fixed in this project
  //   (x: wall-bounded, otherwise periodic)
  const bool periods[NDIMS] = {false, true, true};
  // global array size, real
  const size_t r_gl_sizes[NDIMS] = {glsizes[0], glsizes[1], glsizes[2]};
  // global array size, complex (after-FFTed value in y)
  const size_t c_gl_sizes[NDIMS] = {glsizes[0], glsizes[1] / 2 + 1, glsizes[2]};
  // number of processes in each dimension,
  //   which is to be optimised
  size_t dims_optimum[NDIMS] = {0, 0, 0};
  double wtime_optimum = DBL_MAX;
  const int root = 0;
  int nprocs = 0;
  int myrank = root;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  const size_t nynz = (size_t)nprocs;
  // factorise and decide dims
  for(size_t ny = 1; ny <= nynz; ny++){
    if(0 != nynz % ny){
      // prime decomposition failed
      continue;
    }
    const size_t nz = nynz / ny;
    const size_t dims[NDIMS] = {1, ny, nz};
    // sanitise, for all dimensions, dims should not
    //   exceed the number of grid points
    // NOTE: refer to the complex array,
    //   which is smaller in general
    bool valid = true;
    for(size_t dim1 = 0; dim1 < NDIMS; dim1++){
      const size_t glsize = c_gl_sizes[dim1];
      for(size_t dim0 = 0; dim0 < NDIMS; dim0++){
        const size_t np = dims[dim0];
        if(np > glsize){
          // number of processes is
          //   greater than number of grids
          valid = false;
        }
      }
    }
    if(!valid){
      continue;
    }
    // execute transposes which are used to solve Poisson equation
    //   and check how long they take in total
    // for the time being only transposes when solving Poisson equation are considered,
    //   i.e. implicity and implicitz also request transposes, which are neglected
    sdecomp_info_t * info = NULL;
    if(0 != sdecomp.construct(
          MPI_COMM_WORLD,
          NDIMS,
          (size_t [NDIMS]){dims[0], dims[1], dims[2]},
          periods,
          &info
    )) return 1;
    // initialise pencils and rotations
    double       * r_x1pcnl = NULL;
    double       * r_y1pcnl = NULL;
    fftw_complex * c_y1pcnl = NULL;
    fftw_complex * c_z1pcnl = NULL;
    fftw_complex * c_x2pcnl = NULL;
    sdecomp_transpose_plan_t * r_x1_to_y1 = NULL;
    sdecomp_transpose_plan_t * r_y1_to_x1 = NULL;
    sdecomp_transpose_plan_t * c_y1_to_z1 = NULL;
    sdecomp_transpose_plan_t * c_z1_to_y1 = NULL;
    sdecomp_transpose_plan_t * c_z1_to_x2 = NULL;
    sdecomp_transpose_plan_t * c_x2_to_z1 = NULL;
    if(0 != sdecomp.transpose.construct(info, SDECOMP_X1PENCIL, SDECOMP_Y1PENCIL, r_gl_sizes, sizeof(      double), &r_x1_to_y1)) return 1;
    if(0 != sdecomp.transpose.construct(info, SDECOMP_Y1PENCIL, SDECOMP_X1PENCIL, r_gl_sizes, sizeof(      double), &r_y1_to_x1)) return 1;
    if(0 != sdecomp.transpose.construct(info, SDECOMP_Y1PENCIL, SDECOMP_Z1PENCIL, c_gl_sizes, sizeof(fftw_complex), &c_y1_to_z1)) return 1;
    if(0 != sdecomp.transpose.construct(info, SDECOMP_Z1PENCIL, SDECOMP_Y1PENCIL, c_gl_sizes, sizeof(fftw_complex), &c_z1_to_y1)) return 1;
    if(0 != sdecomp.transpose.construct(info, SDECOMP_Z1PENCIL, SDECOMP_X2PENCIL, c_gl_sizes, sizeof(fftw_complex), &c_z1_to_x2)) return 1;
    if(0 != sdecomp.transpose.construct(info, SDECOMP_X2PENCIL, SDECOMP_Z1PENCIL, c_gl_sizes, sizeof(fftw_complex), &c_x2_to_z1)) return 1;
    size_t r_x1sizes[NDIMS] = {0};
    size_t r_y1sizes[NDIMS] = {0};
    size_t c_y1sizes[NDIMS] = {0};
    size_t c_z1sizes[NDIMS] = {0};
    size_t c_x2sizes[NDIMS] = {0};
    for(sdecomp_dir_t dim = 0; dim < NDIMS; dim++){
      if(0 != sdecomp.get_pencil_mysize(info, SDECOMP_X1PENCIL, dim, r_gl_sizes[dim], r_x1sizes + dim)) return 1;
      if(0 != sdecomp.get_pencil_mysize(info, SDECOMP_Y1PENCIL, dim, r_gl_sizes[dim], r_y1sizes + dim)) return 1;
      if(0 != sdecomp.get_pencil_mysize(info, SDECOMP_Y1PENCIL, dim, c_gl_sizes[dim], c_y1sizes + dim)) return 1;
      if(0 != sdecomp.get_pencil_mysize(info, SDECOMP_Z1PENCIL, dim, c_gl_sizes[dim], c_z1sizes + dim)) return 1;
      if(0 != sdecomp.get_pencil_mysize(info, SDECOMP_X2PENCIL, dim, c_gl_sizes[dim], c_x2sizes + dim)) return 1;
    }
    r_x1pcnl = memory_calloc(r_x1sizes[0] * r_x1sizes[1] * r_x1sizes[2], sizeof(      double));
    r_y1pcnl = memory_calloc(r_y1sizes[0] * r_y1sizes[1] * r_y1sizes[2], sizeof(      double));
    c_y1pcnl = memory_calloc(c_y1sizes[0] * c_y1sizes[1] * c_y1sizes[2], sizeof(fftw_complex));
    c_z1pcnl = memory_calloc(c_z1sizes[0] * c_z1sizes[1] * c_z1sizes[2], sizeof(fftw_complex));
    c_x2pcnl = memory_calloc(c_x2sizes[0] * c_x2sizes[1] * c_x2sizes[2], sizeof(fftw_complex));
    // execute transpose, repeat for "niter" times
    const size_t niter = 4;
    const double tic = timer();
    for(size_t iter = 0; iter < niter; iter++){
      sdecomp.transpose.execute(r_x1_to_y1, r_x1pcnl, r_y1pcnl);
      sdecomp.transpose.execute(c_y1_to_z1, c_y1pcnl, c_z1pcnl);
      sdecomp.transpose.execute(c_z1_to_x2, c_z1pcnl, c_x2pcnl);
      sdecomp.transpose.execute(c_x2_to_z1, c_x2pcnl, c_z1pcnl);
      sdecomp.transpose.execute(c_z1_to_y1, c_z1pcnl, c_y1pcnl);
      sdecomp.transpose.execute(r_y1_to_x1, r_y1pcnl, r_x1pcnl);
    }
    const double toc = timer();
    // clean-up tentative transpose plans and buffers
    sdecomp.transpose.destruct(r_x1_to_y1);
    sdecomp.transpose.destruct(c_y1_to_z1);
    sdecomp.transpose.destruct(c_z1_to_x2);
    sdecomp.transpose.destruct(c_x2_to_z1);
    sdecomp.transpose.destruct(c_z1_to_y1);
    sdecomp.transpose.destruct(r_y1_to_x1);
    memory_free(r_x1pcnl);
    memory_free(r_y1pcnl);
    memory_free(c_y1pcnl);
    memory_free(c_z1pcnl);
    memory_free(c_x2pcnl);
    // clean-up current sdecomp config
    sdecomp.destruct(info);
    // check time
    const double wtime = (toc - tic) / niter;
    if(wtime < wtime_optimum){
      // this is the best option for now,
      //   update candidate
      for(size_t dim = 0; dim < NDIMS; dim++){
        dims_optimum[dim] = dims[dim];
      }
      wtime_optimum = wtime;
    }
    if(root == myrank){
      const int nd = get_ndigits(nprocs);
      printf("dims: [%*zu, %*zu, %*zu]: % .7e [sec]\n", nd, dims[0], nd, dims[1], nd, dims[2], wtime);
    }
  }
  // create sdecomp which will be used in the main run
  if(0 != sdecomp.construct(
        MPI_COMM_WORLD,
        NDIMS,
        (size_t [NDIMS]){dims_optimum[0], dims_optimum[1], dims_optimum[2]},
        periods,
        info_optimum
  )) return 1;
  if(root == myrank){
    const int nd = get_ndigits(nprocs);
    printf("Conclusive domain decomposition: [%*zu, %*zu, %*zu]\n", nd, dims_optimum[0], nd, dims_optimum[1], nd, dims_optimum[2]);
  }
  return 0;
}
#endif

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
  // dxf: distance from cell face to cell face | 8
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
  // dxc: distance from cell center to cell center (generally) | 8
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
#if NDIMS == 3
  double * restrict   dz    = &domain->dz;
#endif
  // load spatial information | 3
  if(0 != domain_load(dirname_ic, domain)){
    return 1;
  }
  // compute grid sizes | 8
  // allocate and initialise x coordinates
  *dxf = allocate_and_init_dxf(glsizes[0], *xf);
  *dxc = allocate_and_init_dxc(glsizes[0], *xc);
  // grid sizes in homogeneous directions
  *dy = lengths[1] / glsizes[1];
#if NDIMS == 3
  *dz = lengths[2] / glsizes[2];
#endif
  // initialise sdecomp to distribute the domain | 14
#if NDIMS == 2
  if(0 != sdecomp.construct(
        MPI_COMM_WORLD,
        NDIMS,
        (size_t [NDIMS]){0, 0},
        (bool [NDIMS]){false, true},
        info
  )) return 1;
#else
  if(0 != optimise_sdecomp_init(
        glsizes,
        info
  )) return 1;
#endif
  // local array sizes and offsets | 4
  for(size_t dim = 0; dim < NDIMS; dim++){
    sdecomp.get_pencil_mysize(*info, SDECOMP_X1PENCIL, dim, glsizes[dim], mysizes + dim);
    sdecomp.get_pencil_offset(*info, SDECOMP_X1PENCIL, dim, glsizes[dim], offsets + dim);
  }
  report(domain);
  return 0;
}

