#if !defined(FLUID_H)
#define FLUID_H

#include "array.h"
#include "domain.h"

// definition of a structure fluid_t_
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
  array_t p;
  array_t psi;
  array_t t;
  array_t srcux[3];
  array_t srcuy[3];
  array_t srct[3];
  double Ra, Pr;
  double m_dif, t_dif;
} fluid_t;

#endif // FLUID_H
