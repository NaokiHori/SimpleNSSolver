#include <math.h>
#include <mpi.h>
#include "domain.h"
#include "fluid.h"
#include "array_macros/domain/dxf.h"
#include "array_macros/domain/dxc.h"
#include "array_macros/fluid/t.h"
#include "internal.h"

static double get_x_contribution(
    const domain_t * domain,
    const fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict dxf = domain->dxf;
  const double * restrict dxc = domain->dxc;
  const double dy = domain->dy;
  const double * restrict t = fluid->t.data;
  const double diffusivity = fluid->t_dif;
  double dissipation = 0.;
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      // x contribution
      const double cellsize = DXF(i  ) * dy;
      // negative direction
      {
        const double dtdx = 1. / DXC(i  ) * (T(i  , j  ) - T(i-1, j  ));
        const double w0 = DXC(i  ) / DXF(i  );
        const double w1 = 1 == i ? 2. : 1.;
        dissipation += diffusivity * 0.5 * w0 * w1 * pow(dtdx, 2.) * cellsize;
      }
      // positive direction
      {
        const double dtdx = 1. / DXC(i+1) * (T(i+1, j  ) - T(i  , j  ));
        const double w0 = DXC(i+1) / DXF(i  );
        const double w1 = isize == i ? 2. : 1.;
        dissipation += diffusivity * 0.5 * w0 * w1 * pow(dtdx, 2.) * cellsize;
      }
    }
  }
  return dissipation;
}

static double get_y_contribution(
    const domain_t * domain,
    const fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict dxf = domain->dxf;
  const double dy = domain->dy;
  const double * restrict t = fluid->t.data;
  const double diffusivity = fluid->t_dif;
  double dissipation = 0.;
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      // y contribution
      const double cellsize = DXF(i  ) * dy;
      // negative direction
      {
        const double dtdy = 1. / dy * (T(i  , j  ) - T(i  , j-1));
        dissipation += diffusivity * 0.5 * pow(dtdy, 2.) * cellsize;
      }
      // positive direction
      {
        const double dtdy = 1. / dy * (T(i  , j+1) - T(i  , j  ));
        dissipation += diffusivity * 0.5 * pow(dtdy, 2.) * cellsize;
      }
    }
  }
  return dissipation;
}

/**
 * @brief compute Nusselt number based on thermal dissipation
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] fluid  : temperature and its diffusivity
 * @return           : Nusselt number
 */
double logging_internal_compute_nu_thermal_energy_dissipation(
    const domain_t * domain,
    const fluid_t * fluid
){
  const int root = 0;
  int myrank = root;
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_rank(domain->info, &myrank);
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  // compute thermal energy dissipation
  // reference heat flux
  const double ref = logging_internal_compute_reference_heat_flux(domain, fluid);
  double retval = 0.;
  retval += get_x_contribution(domain, fluid);
  retval += get_y_contribution(domain, fluid);
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : &retval;
  void * recvbuf = &retval;
  MPI_Reduce(sendbuf, recvbuf, 1, MPI_DOUBLE, MPI_SUM, root, comm_cart);
  retval /= ref;
  return retval;
}

