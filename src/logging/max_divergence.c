#include <math.h>
#include "domain.h"
#include "fluid.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#if NDIMS == 3
#include "array_macros/fluid/uz.h"
#endif
#include "internal.h"

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

/**
 * @brief check maximum divergence and write it to a file
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : domain information
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity
 * @return           : error code
 */
int logging_check_max_divergence (
    const char fname[],
    const domain_t * domain,
    const double time,
    const fluid_t * fluid
) {
  const int root = 0;
  int myrank = root;
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_rank(domain->info, &myrank);
  sdecomp.get_comm_cart(domain->info, &comm_cart);
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
  // check max local divergence | 30
  double divmax = 0.;
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
    // check maximum
    divmax = fmax(divmax, fabs(div));
  END
  // collect information among all processes
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : &divmax;
  void * recvbuf = &divmax;
  MPI_Reduce(sendbuf, recvbuf, 1, MPI_DOUBLE, MPI_MAX, root, comm_cart);
  logging_internal_output(
      fname,
      domain,
      time,
      1,
      &divmax
  );
  return 0;
}

