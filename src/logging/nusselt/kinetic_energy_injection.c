#include <mpi.h>
#include "domain.h"
#include "fluid.h"
#include "array_macros/domain/dxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/t.h"
#include "internal.h"

/**
 * @brief compute Nusselt number based on the total energy injection
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] fluid  : wall-normal velocity (ux) and temperature
 * @return           : Nusselt number
 */
double logging_internal_compute_nu_kinetic_energy_injection(
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
  const double * restrict ux = fluid->ux.data;
  const double * restrict t  = fluid->t .data;
  // energy injection | 35
  // reference heat flux
  const double ref = logging_internal_compute_reference_heat_flux(domain, fluid);
  // integral in the whole domain
  double retval = 0.;
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      const double dx = DXC(i  );
      const double vel = UX(i, j);
      const double t_ = + 0.5 * T(i-1, j) + 0.5 * T(i, j);
      const double integrand = vel * t_;
      const double cellsize = dx * dy;
      retval += integrand * cellsize;
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        const double dx = DXC(i  );
        const double vel = UX(i, j, k);
        const double t_ = + 0.5 * T(i-1, j, k) + 0.5 * T(i, j, k);
        const double integrand = vel * t_;
        const double cellsize = dx * dy * dz;
        retval += integrand * cellsize;
      }
    }
  }
#endif
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : &retval;
  void * recvbuf = &retval;
  MPI_Reduce(sendbuf, recvbuf, 1, MPI_DOUBLE, MPI_SUM, root, comm_cart);
  // normalise by the reference value
  retval /= ref;
  // add laminar contribution
  retval += 1.;
  return retval;
}

