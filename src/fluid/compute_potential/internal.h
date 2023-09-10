#if !defined(FLUID_COMPUTE_POTENTIAL_INTERNAL_H)
#define FLUID_COMPUTE_POTENTIAL_INTERNAL_H

#include "domain.h"
#include "fluid.h"

// compute scalar potential by solving Poisson equation
// (general-purpose version)
extern int fluid_compute_potential_dft(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

// compute scalar potential by solving Poisson equation
// (efficient version)
extern int fluid_compute_potential_dct(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

#endif // FLUID_COMPUTE_POTENTIAL_INTERNAL_H
