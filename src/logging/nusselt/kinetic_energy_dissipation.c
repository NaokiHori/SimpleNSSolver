#include <math.h>
#include <mpi.h>
#include "domain.h"
#include "fluid.h"
#include "array_macros/domain/dxf.h"
#include "array_macros/domain/dxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "internal.h"

static double get_ux_x_contribution(
    const domain_t * domain,
    const fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict dxf = domain->dxf;
  const double * restrict dxc = domain->dxc;
  const double dy = domain->dy;
  const double * restrict ux = fluid->ux.data;
  const double diffusivity = fluid->m_dif;
  double dissipation = 0.;
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize + 1; i++){
      // ux-x contribution
      const double cellsize = DXC(i  ) * dy;
      // negative direction
      if(1 != i){
        const double duxdx = 1. / DXF(i-1) * (UX(i  , j  ) - UX(i-1, j  ));
        const double w = DXF(i-1) / DXC(i  );
        dissipation += diffusivity * 0.5 * w * pow(duxdx, 2.) * cellsize;
      }
      // positive direction
      if(isize + 1 != i){
        const double duxdx = 1. / DXF(i  ) * (UX(i+1, j  ) - UX(i  , j  ));
        const double w = DXF(i  ) / DXC(i  );
        dissipation += diffusivity * 0.5 * w * pow(duxdx, 2.) * cellsize;
      }
    }
  }
  return dissipation;
}

static double get_ux_y_contribution(
    const domain_t * domain,
    const fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict dxc = domain->dxc;
  const double dy = domain->dy;
  const double * restrict ux = fluid->ux.data;
  const double diffusivity = fluid->m_dif;
  double dissipation = 0.;
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize + 1; i++){
      // ux-y contribution
      const double cellsize = DXC(i  ) * dy;
      // negative direction
      {
        const double duxdy = 1. / dy * (UX(i  , j  ) - UX(i  , j-1));
        dissipation += diffusivity * 0.5 * pow(duxdy, 2.) * cellsize;
      }
      // positive direction
      {
        const double duxdy = 1. / dy * (UX(i  , j+1) - UX(i  , j  ));
        dissipation += diffusivity * 0.5 * pow(duxdy, 2.) * cellsize;
      }
    }
  }
  return dissipation;
}

static double get_uy_x_contribution(
    const domain_t * domain,
    const fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict dxf = domain->dxf;
  const double * restrict dxc = domain->dxc;
  const double dy = domain->dy;
  const double * restrict uy = fluid->uy.data;
  const double diffusivity = fluid->m_dif;
  double dissipation = 0.;
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      // uy-x contribution
      const double cellsize = DXF(i  ) * dy;
      // negative direction
      {
        const double duydx = 1. / DXC(i  ) * (UY(i  , j  ) - UY(i-1, j  ));
        const double w0 = DXC(i  ) / DXF(i  );
        const double w1 = 1 == i ? 2. : 1.;
        dissipation += diffusivity * 0.5 * w0 * w1 * pow(duydx, 2.) * cellsize;
      }
      // positive direction
      {
        const double duydx = 1. / DXC(i+1) * (UY(i+1, j  ) - UY(i  , j  ));
        const double w0 = DXC(i+1) / DXF(i  );
        const double w1 = isize == i ? 2. : 1.;
        dissipation += diffusivity * 0.5 * w0 * w1 * pow(duydx, 2.) * cellsize;
      }
    }
  }
  return dissipation;
}

static double get_uy_y_contribution(
    const domain_t * domain,
    const fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict dxf = domain->dxf;
  const double dy = domain->dy;
  const double * restrict uy = fluid->uy.data;
  const double diffusivity = fluid->m_dif;
  double dissipation = 0.;
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      // uy-y contribution
      const double cellsize = DXF(i  ) * dy;
      // negative direction
      {
        const double duydy = 1. / dy * (UY(i  , j  ) - UY(i  , j-1));
        dissipation += diffusivity * 0.5 * pow(duydy, 2.) * cellsize;
      }
      // positive direction
      {
        const double duydy = 1. / dy * (UY(i  , j+1) - UY(i  , j  ));
        dissipation += diffusivity * 0.5 * pow(duydy, 2.) * cellsize;
      }
    }
  }
  return dissipation;
}

/**
 * @brief compute Nusselt number based on kinetic dissipation
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] fluid  : flow field
 * @return           : Nusselt number
 */
double logging_internal_compute_nu_kinetic_energy_dissipation(
    const domain_t * domain,
    const fluid_t * fluid
){
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  // compute kinetic energy dissipation
  // reference heat flux
  const double ref = logging_internal_compute_reference_heat_flux(domain, fluid);
  double retval = 0.;
  // ux contribution
  retval += get_ux_x_contribution(domain, fluid);
  retval += get_ux_y_contribution(domain, fluid);
  // uy contribution
  retval += get_uy_x_contribution(domain, fluid);
  retval += get_uy_y_contribution(domain, fluid);
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : &retval;
  void * recvbuf = &retval;
  MPI_Reduce(sendbuf, recvbuf, 1, MPI_DOUBLE, MPI_SUM, root, comm_cart);
  retval /= ref;
  return retval + 1.;
}

