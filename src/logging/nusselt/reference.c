#include "domain.h"
#include "fluid.h"

double logging_internal_compute_reference_heat_flux(
    const domain_t * domain,
    const fluid_t * fluid
){
  // compute laminar heat flux | 10
  const double ly = domain->lengths[1];
#if NDIMS == 3
  const double lz = domain->lengths[2];
#endif
  const double diffusivity = fluid->t_dif;
#if NDIMS == 2
  return diffusivity * ly;
#else
  return diffusivity * ly * lz;
#endif
}

