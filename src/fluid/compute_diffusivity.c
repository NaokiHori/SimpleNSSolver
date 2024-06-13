#include <math.h>
#include "fluid.h"

double fluid_compute_momentum_diffusivity (
    const fluid_t * fluid
) {
  return sqrt(fluid->Pr) / sqrt(fluid->Ra);
}

double fluid_compute_temperature_diffusivity (
    const fluid_t * fluid
) {
  return 1. / sqrt(fluid->Pr) / sqrt(fluid->Ra);
}

