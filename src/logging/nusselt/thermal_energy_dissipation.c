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
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict dxf = domain->dxf;
  const double * restrict dxc = domain->dxc;
  const double dy = domain->dy;
#if NDIMS == 3
  const double dz = domain->dz;
#endif
  const double * restrict t = fluid->t.data;
  const double diffusivity = fluid->t_dif;
  double dissipation = 0.;
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      // x contribution | 15
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
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // x contribution | 15
        const double cellsize = DXF(i  ) * dy * dz;
        // negative direction
        {
          const double dtdx = 1. / DXC(i  ) * (T(i  , j  , k  ) - T(i-1, j  , k  ));
          const double w0 = DXC(i  ) / DXF(i  );
          const double w1 = 1 == i ? 2. : 1.;
          dissipation += diffusivity * 0.5 * w0 * w1 * pow(dtdx, 2.) * cellsize;
        }
        // positive direction
        {
          const double dtdx = 1. / DXC(i+1) * (T(i+1, j  , k  ) - T(i  , j  , k  ));
          const double w0 = DXC(i+1) / DXF(i  );
          const double w1 = isize == i ? 2. : 1.;
          dissipation += diffusivity * 0.5 * w0 * w1 * pow(dtdx, 2.) * cellsize;
        }
      }
    }
  }
#endif
  return dissipation;
}

static double get_y_contribution(
    const domain_t * domain,
    const fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict dxf = domain->dxf;
  const double dy = domain->dy;
#if NDIMS == 3
  const double dz = domain->dz;
#endif
  const double * restrict t = fluid->t.data;
  const double diffusivity = fluid->t_dif;
  double dissipation = 0.;
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      // y contribution | 11
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
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // y contribution | 11
        const double cellsize = DXF(i  ) * dy * dz;
        // negative direction
        {
          const double dtdy = 1. / dy * (T(i  , j  , k  ) - T(i  , j-1, k  ));
          dissipation += diffusivity * 0.5 * pow(dtdy, 2.) * cellsize;
        }
        // positive direction
        {
          const double dtdy = 1. / dy * (T(i  , j+1, k  ) - T(i  , j  , k  ));
          dissipation += diffusivity * 0.5 * pow(dtdy, 2.) * cellsize;
        }
      }
    }
  }
#endif
  return dissipation;
}

#if NDIMS == 3
static double get_z_contribution(
    const domain_t * domain,
    const fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict dxf = domain->dxf;
  const double dy = domain->dy;
  const double dz = domain->dz;
  const double * restrict t = fluid->t.data;
  const double diffusivity = fluid->t_dif;
  double dissipation = 0.;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // z contribution | 11
        const double cellsize = DXF(i  ) * dy * dz;
        // negative direction
        {
          const double dtdz = 1. / dz * (T(i  , j  , k  ) - T(i  , j  , k-1));
          dissipation += diffusivity * 0.5 * pow(dtdz, 2.) * cellsize;
        }
        // positive direction
        {
          const double dtdz = 1. / dz * (T(i  , j  , k+1) - T(i  , j  , k  ));
          dissipation += diffusivity * 0.5 * pow(dtdz, 2.) * cellsize;
        }
      }
    }
  }
  return dissipation;
}
#endif

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
  // compute thermal energy dissipation | 11
  // reference heat flux
  const double ref = logging_internal_compute_reference_heat_flux(domain, fluid);
  double retval = 0.;
  retval += get_x_contribution(domain, fluid);
  retval += get_y_contribution(domain, fluid);
#if NDIMS == 3
  retval += get_z_contribution(domain, fluid);
#endif
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : &retval;
  void * recvbuf = &retval;
  MPI_Reduce(sendbuf, recvbuf, 1, MPI_DOUBLE, MPI_SUM, root, comm_cart);
  retval /= ref;
  return retval;
}

