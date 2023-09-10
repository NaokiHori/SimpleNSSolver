#include <mpi.h>
#include "domain.h"
#include "fluid.h"
#include "array_macros/domain/dxc.h"
#include "array_macros/fluid/t.h"
#include "internal.h"

/**
 * @brief compute Nusselt number based on the heat flux on the walls
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] fluid  : temperature and diffusivity
 * @return           : Nusselt number
 */
double logging_internal_compute_nu_heat_flux(
    const domain_t * domain,
    const fluid_t * fluid
){
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict dxc = domain->dxc;
  const double dy = domain->dy;
#if NDIMS == 3
  const double dz = domain->dz;
#endif
  const double * restrict t = fluid->t.data;
  const double diffusivity = fluid->t_dif;
  // heat flux on the walls | 33
  // reference heat flux
  const double ref = logging_internal_compute_reference_heat_flux(domain, fluid);
  // integral on the two walls
  double retval = 0.;
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    const double ds = dy;
    const double dx_xm = DXC(      1);
    const double dx_xp = DXC(isize+1);
    const double dt_xm = + T(      1, j) - T(    0, j);
    const double dt_xp = + T(isize+1, j) - T(isize, j);
    const double dtdx_xm = dt_xm / dx_xm;
    const double dtdx_xp = dt_xp / dx_xp;
    // average two walls
    retval -= 0.5 * diffusivity * (dtdx_xm + dtdx_xp) * ds;
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      const double ds = dy * dz;
      const double dx_xm = DXC(      1);
      const double dx_xp = DXC(isize+1);
      const double dt_xm = + T(      1, j, k) - T(    0, j, k);
      const double dt_xp = + T(isize+1, j, k) - T(isize, j, k);
      const double dtdx_xm = dt_xm / dx_xm;
      const double dtdx_xp = dt_xp / dx_xp;
      // average two walls
      retval -= 0.5 * diffusivity * (dtdx_xm + dtdx_xp) * ds;
    }
  }
#endif
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : &retval;
  void * recvbuf = &retval;
  MPI_Reduce(sendbuf, recvbuf, 1, MPI_DOUBLE, MPI_SUM, root, comm_cart);
  return retval / ref;
}

