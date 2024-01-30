#include "runge_kutta.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "internal.h"

/**
 * @brief correct non-solenoidal velocity using scalar potential psi
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : scalar potential (in), velocity (out)
 * @return               : error code
 */
int fluid_correct_velocity(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
){
  // compute prefactor gamma dt
  const double gamma = rkcoefs[rkstep][rk_g];
  const double prefactor = gamma * dt;
  if(0 != fluid_correct_velocity_ux(domain, prefactor, fluid)) return 1;
  if(0 != fluid_correct_velocity_uy(domain, prefactor, fluid)) return 1;
  return 0;
}

