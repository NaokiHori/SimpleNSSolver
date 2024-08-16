#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "sdecomp.h"
#include "memory.h"
#include "runge_kutta.h"
#include "domain.h"
#include "tdm.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "array_macros/domain/xc.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#if NDIMS == 3
#include "array_macros/fluid/uz.h"
#endif
#include "array_macros/fluid/psi.h"

static const double g_pi = 3.14159265358979324;

// structure only used to solve Poisson equation
typedef struct {
  bool is_initialised;
  void * restrict buf0;
  void * restrict buf1;
  fftw_plan fftw_plan_y[2];
#if NDIMS == 3
  fftw_plan fftw_plan_z[2];
#endif
  size_t tdm_sizes[NDIMS];
  tdm_info_t * tdm_info;
  double * restrict eval_ys;
#if NDIMS == 3
  double * restrict eval_zs;
#endif
  sdecomp_transpose_plan_t * r_transposer_x1_to_y1;
  sdecomp_transpose_plan_t * r_transposer_y1_to_x1;
#if NDIMS == 2
  sdecomp_transpose_plan_t * c_transposer_x1_to_y1;
  sdecomp_transpose_plan_t * c_transposer_y1_to_x1;
#else
  sdecomp_transpose_plan_t * c_transposer_y1_to_z1;
  sdecomp_transpose_plan_t * c_transposer_z1_to_y1;
  sdecomp_transpose_plan_t * c_transposer_z1_to_x2;
  sdecomp_transpose_plan_t * c_transposer_x2_to_z1;
#endif
} poisson_solver_t;

/* initialise Poisson solver */
// several pencils for different data types are treated
//   and thus this source is very complicated
// used prefixes are as follows:
//   r_: real    (double)       type
//   c_: complex (fftw_complex) type
//
//   gl_    : global array size (not pencils)
//   x1pncl_: each pencl (x1, y1, ...)

/* size of domain and pencils */
// NOTE: define globally to reduce the number of arguments of static functions
// global domain size in real space
static size_t r_gl_sizes[NDIMS] = {0};
// global domain size in complex space
static size_t c_gl_sizes[NDIMS] = {0};
// local domain size (x1 pencil) in real space
static size_t r_x1pncl_sizes[NDIMS] = {0};
// local domain size (y1 pencil) in real space
static size_t r_y1pncl_sizes[NDIMS] = {0};
// local domain size (y1 pencil) in complex space
static size_t c_y1pncl_sizes[NDIMS] = {0};
#if NDIMS == 2
// local domain size (x1 pencil) in complex space
static size_t c_x1pncl_sizes[NDIMS] = {0};
#else
// local domain size (z1 pencil) in complex space
static size_t c_z1pncl_sizes[NDIMS] = {0};
// local domain size (x2 pencil) in complex space
static size_t c_x2pncl_sizes[NDIMS] = {0};
#endif

static size_t prod (
    const size_t sizes[NDIMS]
) {
  // compute the product of the given vector
  size_t nitems = 1;
  for (size_t dim = 0; dim < NDIMS; dim++) {
    nitems *= sizes[dim];
  }
  return nitems;
}

static int report_failure (
    const char type[]
) {
  // function to just dump error message and abort
  FILE * stream = stderr;
  fprintf(stream, "Poisson solver, initialisation failed: %s\n", type);
  fprintf(stream, "  FFTW:    A possible reason is you link Intel-MKL lib\n");
  fprintf(stream, "           Make sure you use FFTW3 directly,\n");
  fprintf(stream, "           NOT its crazy wrapper offered by MKL\n");
  fprintf(stream, "  SDECOMP: Check sdecomp.log and check arguments\n");
  fprintf(stream, "           If they are all correct, PLEASE CONTACT ME\n");
  fflush(stream);
  return 0;
}

static size_t max (
    const size_t val0,
    const size_t val1
) {
  if (val0 > val1) {
    return val0;
  }else{
    return val1;
  }
}

static int compute_pencil_sizes (
    const domain_t * domain
) {
  // NOTE: those variables are defined globally at the top of this file
  //   to reduce the number of arguments which functions take
  // global domain size in real space
  const sdecomp_info_t * info = domain->info;
  r_gl_sizes[0] = domain->glsizes[0];
  r_gl_sizes[1] = domain->glsizes[1];
#if NDIMS == 3
  r_gl_sizes[2] = domain->glsizes[2];
#endif
  // global domain size in complex space
  // NOTE: Hermite symmetry in y
  c_gl_sizes[0] = domain->glsizes[0];
  c_gl_sizes[1] = domain->glsizes[1] / 2 + 1;
#if NDIMS == 3
  c_gl_sizes[2] = domain->glsizes[2];
#endif
  // local domain sizes
  for (sdecomp_dir_t dim = 0; dim < NDIMS; dim++) {
    if (0 != sdecomp.get_pencil_mysize(info, SDECOMP_X1PENCIL, dim, r_gl_sizes[dim], r_x1pncl_sizes + dim)) return 1;
    if (0 != sdecomp.get_pencil_mysize(info, SDECOMP_Y1PENCIL, dim, r_gl_sizes[dim], r_y1pncl_sizes + dim)) return 1;
    if (0 != sdecomp.get_pencil_mysize(info, SDECOMP_Y1PENCIL, dim, c_gl_sizes[dim], c_y1pncl_sizes + dim)) return 1;
#if NDIMS == 2
    if (0 != sdecomp.get_pencil_mysize(info, SDECOMP_X1PENCIL, dim, c_gl_sizes[dim], c_x1pncl_sizes + dim)) return 1;
#else
    if (0 != sdecomp.get_pencil_mysize(info, SDECOMP_Z1PENCIL, dim, c_gl_sizes[dim], c_z1pncl_sizes + dim)) return 1;
    if (0 != sdecomp.get_pencil_mysize(info, SDECOMP_X2PENCIL, dim, c_gl_sizes[dim], c_x2pncl_sizes + dim)) return 1;
#endif
  }
  return 0;
}

static int allocate_buffers (
    poisson_solver_t * poisson_solver
) {
  // although there are bunch of pencils involved,
  //   two buffers are enough to do the job,
  //   which are allocated here
  void * restrict * buf0 = &poisson_solver->buf0;
  void * restrict * buf1 = &poisson_solver->buf1;
  const size_t r_dsize = sizeof(double);
  const size_t c_dsize = sizeof(fftw_complex);
  size_t buf0_bytes = 0;
  size_t buf1_bytes = 0;
#if NDIMS == 2
  // r_x1pncl -> rotate -> r_y1pncl -> FFT -> c_y1pncl -> rotate -> c_x1pncl
  // buffer0               buffer1            buffer0               buffer1
  buf0_bytes = max(buf0_bytes, r_dsize * prod(r_x1pncl_sizes));
  buf0_bytes = max(buf0_bytes, c_dsize * prod(c_y1pncl_sizes));
  buf1_bytes = max(buf1_bytes, r_dsize * prod(r_y1pncl_sizes));
  buf1_bytes = max(buf1_bytes, c_dsize * prod(c_x1pncl_sizes));
#else
  // r_x1pncl -> rotate -> r_y1pncl -> FFT -> c_y1pncl -> rotate -> c_z1pncl -> FFT -> c_z1pncl -> rotate -> c_x2pncl
  // buffer0               buffer1            buffer0               buffer1            buffer0               buffer1
  buf0_bytes = max(buf0_bytes, r_dsize * prod(r_x1pncl_sizes));
  buf0_bytes = max(buf0_bytes, c_dsize * prod(c_y1pncl_sizes));
  buf0_bytes = max(buf0_bytes, c_dsize * prod(c_z1pncl_sizes));
  buf1_bytes = max(buf1_bytes, r_dsize * prod(r_y1pncl_sizes));
  buf1_bytes = max(buf1_bytes, c_dsize * prod(c_z1pncl_sizes));
  buf1_bytes = max(buf1_bytes, c_dsize * prod(c_x2pncl_sizes));
#endif
  // allocate them using fftw_malloc to enforce them 16bit-aligned for SIMD
  *buf0 = fftw_malloc(buf0_bytes);
  if (NULL == *buf0) {
    fprintf(stderr, "FATAL: fftw_malloc failed (requested %zu bytes)\n", buf0_bytes);
    fflush(stderr);
    return 1;
  }
  *buf1 = fftw_malloc(buf1_bytes);
  if (NULL == *buf1) {
    fprintf(stderr, "FATAL: fftw_malloc failed (requested %zu bytes)\n", buf1_bytes);
    fflush(stderr);
    return 1;
  }
  return 0;
}

static int init_tri_diagonal_solver (
    const domain_t * domain,
    poisson_solver_t * poisson_solver
) {
  // tri-diagonal matrice, whose size is Nx,
  //   are solved for Ny x Nz times
  // since lower- and upper-diagonal components are
  //   independent to y, z directions and time,
  //   we compute here and re-use them
  // center-diagonal components are, on the other hand,
  //   dependent on time and thus needs to compute everytime
  //   in the solver
  size_t * restrict tdm_sizes = poisson_solver->tdm_sizes;
  tdm_info_t ** tdm_info = &poisson_solver->tdm_info;
#if NDIMS == 2
  tdm_sizes[0] = c_x1pncl_sizes[0];
  tdm_sizes[1] = c_x1pncl_sizes[1];
#else
  tdm_sizes[0] = c_x2pncl_sizes[0];
  tdm_sizes[1] = c_x2pncl_sizes[1];
  tdm_sizes[2] = c_x2pncl_sizes[2];
#endif
  if (0 != tdm.construct(
    /* size of system */ tdm_sizes[0],
    /* number of rhs  */ 1,
    /* is periodic    */ false,
    /* is complex     */ true,
    /* output         */ tdm_info
  )) return 1;
  // initialise tri-diagonal matrix in x direction | 12
  double * tdm_l = NULL;
  double * tdm_u = NULL;
  tdm.get_l(*tdm_info, &tdm_l);
  tdm.get_u(*tdm_info, &tdm_u);
  const double * hxxf = domain->hxxf;
  const double * jdxf = domain->jdxf;
  const double * jdxc = domain->jdxc;
  for (int i = 1; i <= (int)tdm_sizes[0]; i++) {
    // N.B. iterated from i = 1 to use macros
    tdm_l[i-1] = 1. / JDXC(i  ) * JDXF(i  ) / HXXF(i  ) / HXXF(i  );
    tdm_u[i-1] = 1. / JDXC(i  ) * JDXF(i+1) / HXXF(i+1) / HXXF(i+1);
  }
  return 0;
}

static int init_pencil_rotations (
    const domain_t * domain,
    poisson_solver_t * poisson_solver
) {
  const sdecomp_info_t * info = domain->info;
  const size_t r_dsize = sizeof(double);
  const size_t c_dsize = sizeof(fftw_complex);
  if (0 != sdecomp.transpose.construct(info, SDECOMP_X1PENCIL, SDECOMP_Y1PENCIL, r_gl_sizes, r_dsize, &poisson_solver->r_transposer_x1_to_y1)) {
    report_failure("SDECOMP x1 to y1 for real");
    return 1;
  }
  if (0 != sdecomp.transpose.construct(info, SDECOMP_Y1PENCIL, SDECOMP_X1PENCIL, r_gl_sizes, r_dsize, &poisson_solver->r_transposer_y1_to_x1)) {
    report_failure("SDECOMP y1 to x1 for real");
    return 1;
  }
#if NDIMS == 2
  if (0 != sdecomp.transpose.construct(info, SDECOMP_Y1PENCIL, SDECOMP_X1PENCIL, c_gl_sizes, c_dsize, &poisson_solver->c_transposer_y1_to_x1)) {
    report_failure("SDECOMP y1 to x1 for complex");
    return 1;
  }
  if (0 != sdecomp.transpose.construct(info, SDECOMP_X1PENCIL, SDECOMP_Y1PENCIL, c_gl_sizes, c_dsize, &poisson_solver->c_transposer_x1_to_y1)) {
    report_failure("SDECOMP x1 to y1 for complex");
    return 1;
  }
#else
  if (0 != sdecomp.transpose.construct(info, SDECOMP_Y1PENCIL, SDECOMP_Z1PENCIL, c_gl_sizes, c_dsize, &poisson_solver->c_transposer_y1_to_z1)) {
    report_failure("SDECOMP y1 to z1 for complex");
    return 1;
  }
  if (0 != sdecomp.transpose.construct(info, SDECOMP_Z1PENCIL, SDECOMP_Y1PENCIL, c_gl_sizes, c_dsize, &poisson_solver->c_transposer_z1_to_y1)) {
    report_failure("SDECOMP z1 to y1 for complex");
    return 1;
  }
  if (0 != sdecomp.transpose.construct(info, SDECOMP_Z1PENCIL, SDECOMP_X2PENCIL, c_gl_sizes, c_dsize, &poisson_solver->c_transposer_z1_to_x2)) {
    report_failure("SDECOMP z1 to x2 for complex");
    return 1;
  }
  if (0 != sdecomp.transpose.construct(info, SDECOMP_X2PENCIL, SDECOMP_Z1PENCIL, c_gl_sizes, c_dsize, &poisson_solver->c_transposer_x2_to_z1)) {
    report_failure("SDECOMP x2 to z1 for complex");
    return 1;
  }
#endif
  return 0;
}

static int init_ffts (
    poisson_solver_t * poisson_solver
) {
  const unsigned flags = FFTW_MEASURE;
  // NOTE: two buffers should be properly given
  //   see "allocate_buffers" above
  // y, real / complex
  {
    fftw_plan * fplan = &poisson_solver->fftw_plan_y[0];
    fftw_plan * bplan = &poisson_solver->fftw_plan_y[1];
    const int r_signal_length = r_y1pncl_sizes[SDECOMP_YDIR];
    const int c_signal_length = c_y1pncl_sizes[SDECOMP_YDIR];
#if NDIMS == 2
    const int repeat_for = r_y1pncl_sizes[SDECOMP_XDIR];
#else
    const int repeat_for = r_y1pncl_sizes[SDECOMP_ZDIR] * r_y1pncl_sizes[SDECOMP_XDIR];
#endif
    *fplan = fftw_plan_many_dft_r2c(
        1, &r_signal_length, repeat_for,
        poisson_solver->buf1, NULL, 1, r_signal_length,
        poisson_solver->buf0, NULL, 1, c_signal_length,
        flags
    );
    *bplan = fftw_plan_many_dft_c2r(
        1, &r_signal_length, repeat_for,
        poisson_solver->buf0, NULL, 1, c_signal_length,
        poisson_solver->buf1, NULL, 1, r_signal_length,
        flags
    );
    if (NULL == *fplan) {
      report_failure("FFTW y-forward");
      return 1;
    }
    if (NULL == *bplan) {
      report_failure("FFTW y-backward");
      return 1;
    }
  }
#if NDIMS == 3
  // z, complex to complex
  {
    const int signal_length = c_z1pncl_sizes[SDECOMP_ZDIR];
    const int repeat_for = c_z1pncl_sizes[SDECOMP_XDIR] * c_z1pncl_sizes[SDECOMP_YDIR];
    fftw_plan * fplan = &poisson_solver->fftw_plan_z[0];
    fftw_plan * bplan = &poisson_solver->fftw_plan_z[1];
    *fplan = fftw_plan_many_dft(
        1, &signal_length, repeat_for,
        poisson_solver->buf1, NULL, 1, signal_length,
        poisson_solver->buf0, NULL, 1, signal_length,
        FFTW_FORWARD, flags
    );
    *bplan = fftw_plan_many_dft(
        1, &signal_length, repeat_for,
        poisson_solver->buf0, NULL, 1, signal_length,
        poisson_solver->buf1, NULL, 1, signal_length,
        FFTW_BACKWARD, flags
    );
    if (NULL == *fplan) {
      report_failure("FFTW z-forward");
      return 1;
    }
    if (NULL == *bplan) {
      report_failure("FFTW z-backward");
      return 1;
    }
  }
#endif
  return 0;
}

static int init_eigenvalues (
    const domain_t * domain,
    poisson_solver_t * poisson_solver
) {
  const sdecomp_info_t * info = domain->info;
#if NDIMS == 2
  // x1 pencil
  const sdecomp_pencil_t pencil = SDECOMP_X1PENCIL;
#else
  // x2 pencil
  const sdecomp_pencil_t pencil = SDECOMP_X2PENCIL;
#endif
  // y eigenvalues | 16
  {
    double * restrict * evals = &poisson_solver->eval_ys;
    size_t mysize = 0;
    size_t myoffs = 0;
    sdecomp.get_pencil_mysize(info, pencil, SDECOMP_YDIR, c_gl_sizes[SDECOMP_YDIR], &mysize);
    sdecomp.get_pencil_offset(info, pencil, SDECOMP_YDIR, c_gl_sizes[SDECOMP_YDIR], &myoffs);
    *evals = memory_calloc(mysize, sizeof(double));
    for (size_t local_n = 0; local_n < mysize; local_n++) {
      const size_t global_n = local_n + myoffs;
      (*evals)[local_n] =
        - 4. * pow(
          sin( g_pi * global_n / r_gl_sizes[SDECOMP_YDIR] ),
          2.
        );
    }
  }
#if NDIMS == 3
  // z eigenvalues | 16
  {
    double * restrict * evals = &poisson_solver->eval_zs;
    size_t mysize = 0;
    size_t myoffs = 0;
    sdecomp.get_pencil_mysize(info, pencil, SDECOMP_ZDIR, c_gl_sizes[SDECOMP_ZDIR], &mysize);
    sdecomp.get_pencil_offset(info, pencil, SDECOMP_ZDIR, c_gl_sizes[SDECOMP_ZDIR], &myoffs);
    *evals = memory_calloc(mysize, sizeof(double));
    for (size_t local_n = 0; local_n < mysize; local_n++) {
      const size_t global_n = local_n + myoffs;
      (*evals)[local_n] =
        - 4. * pow(
          sin( g_pi * global_n / r_gl_sizes[SDECOMP_ZDIR] ),
          2.
        );
    }
  }
#endif
  return 0;
}

static int init_poisson_solver (
    const domain_t * domain,
    poisson_solver_t * poisson_solver
) {
  // check domain size (global, local, pencils)
  if (0 != compute_pencil_sizes(domain)) return 1;
  // initialise each part of poisson_solver_t
  if (0 != allocate_buffers(poisson_solver))                 return 1;
  if (0 != init_tri_diagonal_solver(domain, poisson_solver)) return 1;
  if (0 != init_pencil_rotations(domain, poisson_solver))    return 1;
  if (0 != init_ffts(poisson_solver))                        return 1;
  if (0 != init_eigenvalues(domain, poisson_solver))         return 1;
  poisson_solver->is_initialised = true;
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if (root == myrank) {
    printf("DFT-based solver is used\n");
  }
  return 0;
}

#if NDIMS == 2
#define BEGIN \
  for (int cnt = 0, j = 1; j <= jsize; j++) { \
    for (int i = 1; i <= isize; i++, cnt++) {
#define END \
    } \
  }
#else
#define BEGIN \
  for (int cnt = 0, k = 1; k <= ksize; k++) { \
    for (int j = 1; j <= jsize; j++) { \
      for (int i = 1; i <= isize; i++, cnt++) {
#define END \
      } \
    } \
  }
#endif

// compute right-hand side of Poisson equation | 61
static int assign_input (
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    const fluid_t * fluid,
    double * restrict rhs
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxf = domain->hxxf;
  const double hy = domain->hy;
#if NDIMS == 3
  const double hz = domain->hz;
#endif
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict ux = fluid->ux.data;
  const double * restrict uy = fluid->uy.data;
#if NDIMS == 3
  const double * restrict uz = fluid->uz.data;
#endif
  // normalise FFT beforehand
#if NDIMS == 2
  const double norm = 1. * domain->glsizes[1];
#else
  const double norm = 1. * domain->glsizes[1] * domain->glsizes[2];
#endif
  const double prefactor = 1. / (rkcoefs[rkstep].gamma * dt) / norm;
  BEGIN
    const double hx_xm = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i  );
    const double jd_x0 = JDXC(i  );
    const double jd_xp = JDXF(i+1);
#if NDIMS == 2
    const double ux_xm = UX(i  , j  );
    const double ux_xp = UX(i+1, j  );
    const double uy_ym = UY(i  , j  );
    const double uy_yp = UY(i  , j+1);
#else
    const double ux_xm = UX(i  , j  , k  );
    const double ux_xp = UX(i+1, j  , k  );
    const double uy_ym = UY(i  , j  , k  );
    const double uy_yp = UY(i  , j+1, k  );
    const double uz_zm = UZ(i  , j  , k  );
    const double uz_zp = UZ(i  , j  , k+1);
#endif
    const double div = 1. / jd_x0 * (
        - jd_xm / hx_xm * ux_xm + jd_xp / hx_xp * ux_xp
        - jd_x0 / hy    * uy_ym + jd_x0 / hy    * uy_yp
#if NDIMS == 3
        - jd_x0 / hz    * uz_zm + jd_x0 / hz    * uz_zp
#endif
    );
    rhs[cnt] = prefactor * div;
  END
  return 0;
}

static int solve_linear_systems (
    const domain_t * domain,
    poisson_solver_t * poisson_solver
) {
  const double hy = domain->hy;
#if NDIMS == 3
  const double hz = domain->hz;
#endif
  const size_t * restrict tdm_sizes = poisson_solver->tdm_sizes;
  const size_t isize = tdm_sizes[0];
  const size_t jsize = tdm_sizes[1];
#if NDIMS == 3
  const size_t ksize = tdm_sizes[2];
#endif
  // tri-diagonal matrix
  tdm_info_t * tdm_info = poisson_solver->tdm_info;
  double * restrict tdm_l = NULL;
  double * restrict tdm_u = NULL;
  double * restrict tdm_c = NULL;
  tdm.get_l(tdm_info, &tdm_l);
  tdm.get_u(tdm_info, &tdm_u);
  tdm.get_c(tdm_info, &tdm_c);
  // eigenvalues coming from Fourier projection
  const double * restrict eval_ys = poisson_solver->eval_ys;
#if NDIMS == 3
  const double * restrict eval_zs = poisson_solver->eval_zs;
#endif
  fftw_complex * restrict rhs = poisson_solver->buf1;
  // solve tri-diagonal matrices | 37
#if NDIMS == 2
  for (size_t j = 0; j < jsize; j++) {
    const double eval_y = eval_ys[j];
    // set center diagonal components
    for (size_t i = 0; i < isize; i++) {
      tdm_c[i] =
        - tdm_l[i]
        - tdm_u[i]
        + eval_y / hy / hy;
    }
    // boundary treatment (Neumann boundary condition)
    tdm_c[        0] += tdm_l[        0];
    tdm_c[isize - 1] += tdm_u[isize - 1];
    // solve
    tdm.solve(tdm_info, rhs + j * isize);
  }
#else
  for (size_t k = 0; k < ksize; k++) {
    const double eval_z = eval_zs[k];
    for (size_t j = 0; j < jsize; j++) {
      const double eval_y = eval_ys[j];
      // set center diagonal components
      for (size_t i = 0; i < isize; i++) {
        tdm_c[i] =
          - tdm_l[i]
          - tdm_u[i]
          + eval_y / hy / hy
          + eval_z / hz / hz;
      }
      // boundary treatment (Neumann boundary condition)
      tdm_c[        0] += tdm_l[        0];
      tdm_c[isize - 1] += tdm_u[isize - 1];
      // solve
      tdm.solve(tdm_info, rhs + (k * jsize + j) * isize);
    }
  }
#endif
  return 0;
}

static int extract_output (
    const domain_t * domain,
    const double * restrict rhs,
    fluid_t * fluid
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  double * restrict psi = fluid->psi.data;
  BEGIN
#if NDIMS == 2
    PSI(i, j) = rhs[cnt];
#else
    PSI(i, j, k) = rhs[cnt];
#endif
  END
  if (0 != fluid_update_boundaries_psi(domain, &fluid->psi)) {
    return 1;
  }
  return 0;
}

// solve Poisson equation | 105
int fluid_compute_potential (
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
) {
  static poisson_solver_t poisson_solver = {
    .is_initialised = false,
  };
  // initialise Poisson solver
  if (!poisson_solver.is_initialised) {
    if (0 != init_poisson_solver(domain, &poisson_solver)) {
      // failed to initialise Poisson solver
      return 1;
    }
  }
  // compute right-hand side of Poisson equation
  // assigned to buf0
  assign_input(domain, rkstep, dt, fluid, poisson_solver.buf0);
  // transpose real x1pencil to y1pencil
  // from buf0 to buf1
  sdecomp.transpose.execute(
      poisson_solver.r_transposer_x1_to_y1,
      poisson_solver.buf0,
      poisson_solver.buf1
  );
  // project y to wave space
  // f(x, y)    -> f(x, k_y)
  // f(x, y, z) -> f(x, k_y, z)
  // from buf1 to buf0
  fftw_execute(poisson_solver.fftw_plan_y[0]);
#if NDIMS == 2
  // transpose complex y1pencil to x1pencil
  // from buf0 to buf1
  sdecomp.transpose.execute(
      poisson_solver.c_transposer_y1_to_x1,
      poisson_solver.buf0,
      poisson_solver.buf1
  );
#else
  // transpose complex y1pencil to z1pencil
  // from buf0 to buf1
  sdecomp.transpose.execute(
      poisson_solver.c_transposer_y1_to_z1,
      poisson_solver.buf0,
      poisson_solver.buf1
  );
  // project z to wave space
  // f(x, k_y, z) -> f(x, k_y, k_z)
  // from buf1 to buf0
  fftw_execute(poisson_solver.fftw_plan_z[0]);
  // transpose complex z1pencil to x2pencil
  // from buf0 to buf1
  sdecomp.transpose.execute(
      poisson_solver.c_transposer_z1_to_x2,
      poisson_solver.buf0,
      poisson_solver.buf1
  );
#endif
  // solve linear systems in x direction
  solve_linear_systems(domain, &poisson_solver);
#if NDIMS == 2
  // transpose complex x1pencil to y1pencil
  // from buf1 to buf0
  sdecomp.transpose.execute(
      poisson_solver.c_transposer_x1_to_y1,
      poisson_solver.buf1,
      poisson_solver.buf0
  );
#else
  // transpose complex x2pencil to z1pencil
  // from buf1 to buf0
  sdecomp.transpose.execute(
      poisson_solver.c_transposer_x2_to_z1,
      poisson_solver.buf1,
      poisson_solver.buf0
  );
  // project z to physical space
  // f(x, k_y, k_z) -> f(x, k_y, z)
  // from buf0 to buf1
  fftw_execute(poisson_solver.fftw_plan_z[1]);
  // transpose complex z1pencil to y1pencil
  // from buf1 to buf0
  sdecomp.transpose.execute(
      poisson_solver.c_transposer_z1_to_y1,
      poisson_solver.buf1,
      poisson_solver.buf0
  );
#endif
  // project y to physical space
  // f(x, k_y)    -> f(x, y)
  // f(x, k_y, z) -> f(x, y, z)
  // from buf0 to buf1
  fftw_execute(poisson_solver.fftw_plan_y[1]);
  // transpose real y1pencil to x1pencil
  // from buf1 to buf0
  sdecomp.transpose.execute(
      poisson_solver.r_transposer_y1_to_x1,
      poisson_solver.buf1,
      poisson_solver.buf0
  );
  // copy result to the original buffer
  extract_output(domain, poisson_solver.buf0, fluid);
  return 0;
}

