#include <stdbool.h>
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "internal.h"

int fluid_compute_potential(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
){
  // check grid in x direction is uniform
  //   if it is, I can use more efficient algorithm
  //   to compute the scalar potential;
  //   otherwise a versatile version is adopted
  static bool is_initialised = false;
  static bool x_grid_is_uniform = false;
  if(!is_initialised){
    if(0 != domain_check_x_grid_is_uniform(domain, &x_grid_is_uniform)){
      return 1;
    }
    is_initialised = true;
  }
  // use one of
  //   DCT version: fluid_compute_potential_dct
  //   DFT version: fluid_compute_potential_dft
  // depending on the x grid configuration
  if(x_grid_is_uniform){
    if(0 != fluid_compute_potential_dct(domain, rkstep, dt, fluid)){
      return 1;
    }
  }else{
    if(0 != fluid_compute_potential_dft(domain, rkstep, dt, fluid)){
      return 1;
    }
  }
  return 0;
}

