#include "domain.h"
#include "fluid.h"

double logging_internal_compute_reference_heat_flux(
    const domain_t * domain,
    const fluid_t * fluid
){
  // compute laminar heat flux
  const double ly = domain->lengths[1];
  const double lz = domain->lengths[2];
  const double diffusivity = fluid->t_dif;
  return diffusivity * ly * lz;
}

