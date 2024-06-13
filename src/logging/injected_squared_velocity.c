#include "domain.h"
#include "fluid.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/t.h"
#include "internal.h"

#if NDIMS == 2
#define BEGIN \
  for (int j = 1; j <= jsize; j++) { \
    for (int i = 2; i <= isize; i++) {
#define END \
    } \
  }
#else
#define BEGIN \
  for (int k = 1; k <= ksize; k++) { \
    for (int j = 1; j <= jsize; j++) { \
      for (int i = 2; i <= isize; i++) {
#define END \
      } \
    } \
  }
#endif

/**
 * @brief compute injected squared velocity
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity and temperature field
 * @return           : error code
 */
int logging_check_injected_squared_velocity (
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
  const double * restrict jdxf = domain->jdxf;
  const double * restrict ux = fluid->ux.data;
  const double * restrict  t = fluid-> t.data;
  // compute injected squared velocity | 19
  double quantity = 0.;
  BEGIN
    const double jd_x0 = JDXF(i);
#if NDIMS == 2
    const double ux_x0 = UX(i, j);
    const double  t_x0 =
      + 0.5 * T(i-1, j  )
      + 0.5 * T(i  , j  );
#else
    const double ux_x0 = UX(i, j, k);
    const double  t_x0 =
      + 0.5 * T(i-1, j  , k  )
      + 0.5 * T(i  , j  , k  );
#endif
    quantity += jd_x0 * ux_x0 * t_x0;
  END
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : &quantity;
  void * recvbuf = &quantity;
  MPI_Reduce(sendbuf, recvbuf, 1, MPI_DOUBLE, MPI_SUM, root, comm_cart);
  logging_internal_output(
      fname,
      domain,
      time,
      1,
      &quantity
  );
  return 0;
}

