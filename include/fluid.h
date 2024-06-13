#if !defined(FLUID_H)
#define FLUID_H

#include "array.h"
#include "domain.h"

// definition of a structure fluid_t_ | 30
/**
 * @struct fluid_t
 * @brief struct storing fluid-related variables
 * @var ux, uy, uz   : velocity in each direction
 * @var p, psi       : pressure, scalar potential
 * @var t            : temperature
 * @var srcux        : Runge-Kutta source terms for ux
 * @var srcuy        : Runge-Kutta source terms for uy
 * @var srcuz        : Runge-Kutta source terms for uz
 * @var Ra, Pr       : non-dimensional parameters
 * @var m_dif, t_dif : momentum / temperature diffusivities
 */
typedef struct {
  array_t ux;
  array_t uy;
#if NDIMS == 3
  array_t uz;
#endif
  array_t p;
  array_t psi;
  array_t t;
  array_t srcux[3];
  array_t srcuy[3];
#if NDIMS == 3
  array_t srcuz[3];
#endif
  array_t srct[3];
  double Ra, Pr;
} fluid_t;

// initialiser of fluid_t
extern int fluid_init(
    const char dirname_ic[],
    const domain_t * domain,
    fluid_t * fluid
);

// save flow field
extern int fluid_save(
    const char dirname[],
    const domain_t * domain,
    const fluid_t * fluid
);

#endif // FLUID_H
