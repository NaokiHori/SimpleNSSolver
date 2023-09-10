#include "runge_kutta.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "integrate.h"

// integrate the equations for one time step
int integrate(
    const domain_t * domain,
    fluid_t * fluid,
    double * dt
){
  // decide time step size
  if(0 != fluid_decide_dt(domain, fluid, dt)){
    return 1;
  }
  // Runge-Kutta iterations
  // max iteration, should be three
  for(size_t rkstep = 0; rkstep < RKSTEPMAX; rkstep++){
    // predict flow field | 14
    // compute right-hand-side terms of RK scheme
    if(0 != fluid_compute_rhs(domain, fluid)){
      return 1;
    }
    // couple external factors, by default buoyancy force
    if(0 != fluid_couple_external_force(domain, fluid)){
      return 1;
    }
    // update flow field
    // NOTE: while the temperature is fully updated here,
    //   the velocity field is still non-solenoidal
    if(0 != fluid_predict_field(domain, rkstep, *dt, fluid)){
      return 1;
    }
    // compute scalar potential | 6
    // now the temperature field has been updated,
    //   while the velocity field is not divergence free
    //   and thus the following correction step is needed
    if(0 != fluid_compute_potential(domain, rkstep, *dt, fluid)){
      return 1;
    }
    // correct velocity field to satisfy mass conservation | 3
    if(0 != fluid_correct_velocity(domain, rkstep, *dt, fluid)){
      return 1;
    }
    // update pressure | 3
    if(0 != fluid_update_pressure(domain, rkstep, *dt, fluid)){
      return 1;
    }
  }
  return 0;
}

