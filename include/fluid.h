#if !defined(FLUID_H)
#define FLUID_H

#include "array.h"
#include "runge_kutta.h"
#include "domain.h"

// definition of a structure fluid_t
/**
 * @struct fluid_t
 * @brief struct storing fluid-related variables
 * @var u[x-z]    : velocity in each direction
 * @var p, psi    : pressure, scalar potential
 * @var t         : temperature
 * @var srcu[x-z] : Runge-Kutta source terms for ux, uy, uz
 * @var srct      : Runge-Kutta source terms for temperature
 * @var Ra, Pr    : non-dimensional parameters
 */
typedef struct {
  array_t ux;
  array_t uy;
  array_t uz;
  array_t p;
  array_t psi;
  array_t t;
  rkbuffers_t srcux;
  rkbuffers_t srcuy;
  rkbuffers_t srcuz;
  rkbuffers_t srct;
  double Ra;
  double Pr;
} fluid_t;

// initialiser of fluid_t
extern int fluid_init (
    const char dirname_ic[],
    const domain_t * const domain,
    fluid_t * const fluid
);

// save flow field
extern int fluid_save (
    const char dirname[],
    const domain_t * const domain,
    const fluid_t * const fluid
);

#endif // FLUID_H
