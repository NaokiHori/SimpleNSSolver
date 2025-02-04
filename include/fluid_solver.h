#if !defined(FLUID_SOLVER_H)
#define FLUID_SOLVER_H

#include "array.h"
#include "domain.h"
#include "fluid.h"

extern double fluid_compute_momentum_diffusivity (
    const fluid_t * fluid
);

extern double fluid_compute_temperature_diffusivity (
    const fluid_t * fluid
);

// predict the new velocity field and update the temperature field
extern int fluid_predict_field(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

// compute scalar potential by solving Poisson equation
extern int fluid_compute_potential(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

// correct velocity field using scalar potential to enforce divergence zero
extern int fluid_correct_velocity(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

extern int fluid_update_pressure(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

extern int fluid_update_boundaries_ux(
    const domain_t * domain,
    array_t * ux
);

extern int fluid_update_boundaries_uy(
    const domain_t * domain,
    array_t * uy
);

extern int fluid_update_boundaries_p(
    const domain_t * domain,
    array_t * p
);

extern int fluid_update_boundaries_psi(
    const domain_t * domain,
    array_t * psi
);

extern int fluid_update_boundaries_t(
    const domain_t * domain,
    array_t * t
);

#endif // FLUID_SOLVER_H
