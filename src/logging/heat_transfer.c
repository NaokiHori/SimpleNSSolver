#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/fluid/t.h"
#include "internal.h"

#if NDIMS == 2
#define BEGIN \
  for (int j = 1; j <= jsize; j++) {
#define END \
  }
#else
#define BEGIN \
  for (int k = 1; k <= ksize; k++) { \
    for (int j = 1; j <= jsize; j++) {
#define END \
    } \
  }
#endif

/**
 * @brief compute net heat transfer on the walls
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] time   : current simulation time
 * @param[in] fluid  : diffusivity and temperature field
 * @return           : error code
 */
int logging_check_heat_transfer (
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
  const double * restrict jdxf = domain->jdxf;
  const double * restrict t = fluid->t.data;
  // compute heat transfer on the walls | 22
  const double diffusivity = fluid_compute_temperature_diffusivity(fluid);
  // on the bottom and top walls
  double energies[2] = {0., 0.};
  BEGIN
    const double hx_xm = HXXF(        1);
    const double hx_xp = HXXF(isize + 1);
    const double jd_xm = JDXF(        1);
    const double jd_xp = JDXF(isize + 1);
#if NDIMS == 2
    const double dt_xm = - T(    0, j) + T(        1, j);
    const double dt_xp = - T(isize, j) + T(isize + 1, j);
#else
    const double dt_xm = - T(    0, j, k) + T(        1, j, k);
    const double dt_xp = - T(isize, j, k) + T(isize + 1, j, k);
#endif
    energies[0] -= diffusivity * jd_xm / hx_xm / hx_xm * dt_xm;
    energies[1] -= diffusivity * jd_xp / hx_xp / hx_xp * dt_xp;
  END
  const size_t nitems = sizeof(energies) / sizeof(energies[0]);
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : energies;
  void * recvbuf = energies;
  MPI_Reduce(sendbuf, recvbuf, nitems, MPI_DOUBLE, MPI_SUM, root, comm_cart);
  logging_internal_output(
      fname,
      domain,
      time,
      nitems,
      energies
  );
  return 0;
}

