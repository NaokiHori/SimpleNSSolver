#include "domain.h"
#include "fluid.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#if NDIMS == 3
#include "array_macros/fluid/uz.h"
#endif
#include "internal.h"

/**
 * @brief compute total momenta
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : information about domain decomposition and size
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity
 * @return           : error code
 */
int logging_check_total_momentum (
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
  const double * restrict jdxc = domain->jdxc;
  const double * restrict ux = fluid->ux.data;
  const double * restrict uy = fluid->uy.data;
#if NDIMS == 3
  const double * restrict uz = fluid->uz.data;
#endif
  double momenta[NDIMS] = {0.};
  // compute total x-momentum | 15
#if NDIMS == 2
  for (int j = 1; j <= jsize; j++) {
    for (int i = 2; i <= isize; i++) {
      momenta[0] += JDXF(i) * UX(i, j);
    }
  }
#else
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 2; i <= isize; i++) {
        momenta[0] += JDXF(i) * UX(i, j, k);
      }
    }
  }
#endif
  // compute total y-momentum | 15
#if NDIMS == 2
  for (int j = 1; j <= jsize; j++) {
    for (int i = 1; i <= isize; i++) {
      momenta[1] += JDXC(i) * UY(i, j);
    }
  }
#else
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize; i++) {
        momenta[1] += JDXC(i) * UY(i, j, k);
      }
    }
  }
#endif
#if NDIMS == 3
  // compute total z-momentum | 7
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize; i++) {
        momenta[2] += JDXC(i) * UZ(i, j, k);
      }
    }
  }
#endif
  const size_t nitems = sizeof(momenta) / sizeof(momenta[0]);
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : momenta;
  void * recvbuf = momenta;
  MPI_Reduce(sendbuf, recvbuf, nitems, MPI_DOUBLE, MPI_SUM, root, comm_cart);
  logging_internal_output(
      fname,
      domain,
      time,
      nitems,
      momenta
  );
  return 0;
}

