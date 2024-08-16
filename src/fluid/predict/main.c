#include "runge_kutta.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "internal.h"

/**
 * @brief update fields using the previously-computed RK source terms
 * @param[in]     domain : information related to MPI domain decomposition
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : RK source terms (in), flow field (out)
 * @return               : error code
 */
int fluid_predict_field(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
){
  // reset buffers
  // copy previous k-step source term and reset | 14
  if (0 != rkbuffers_reset(&fluid->srcux)) {
    return 1;
  }
  if (0 != rkbuffers_reset(&fluid->srcuy)) {
    return 1;
  }
#if NDIMS == 3
  if (0 != rkbuffers_reset(&fluid->srcuz)) {
    return 1;
  }
#endif
  if (0 != rkbuffers_reset(&fluid->srct)) {
    return 1;
  }
  // compute right-hand-side terms and store them to the corresponding buffers
  if(0 != compute_rhs_ux(domain, fluid)){
    return 1;
  }
  if(0 != compute_rhs_uy(domain, fluid)){
    return 1;
  }
#if NDIMS == 3
  if(0 != compute_rhs_uz(domain, fluid)){
    return 1;
  }
#endif
  if(0 != compute_rhs_t (domain, fluid)){
    return 1;
  }
  // update flow field
  if(0 != update_ux(domain, rkstep, dt, fluid)){
    return 1;
  }
  if(0 != update_uy(domain, rkstep, dt, fluid)){
    return 1;
  }
#if NDIMS == 3
  if(0 != update_uz(domain, rkstep, dt, fluid)){
    return 1;
  }
#endif
  if(0 != update_t (domain, rkstep, dt, fluid)){
    return 1;
  }
  return 0;
}

