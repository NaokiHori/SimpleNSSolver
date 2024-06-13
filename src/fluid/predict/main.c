#include <string.h>
#include "runge_kutta.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "internal.h"

static int reset_srcs(
    array_t * restrict srca,
    array_t * restrict srcb,
    array_t * restrict srcg
){
  // stash previous RK source term,
  //   which is achieved by swapping
  //   the pointers to "data"
  double * tmp = srca->data;
  srca->data = srcb->data;
  srcb->data = tmp;
  // zero-clear current RK source terms (exp/imp)
  // NOTE: when 0 == rkstep, rkcoef for "beta" is zero
  //   and thus zero-clearing"beta" buffers is not needed
  memset(srca->data, 0, srca->datasize);
  memset(srcg->data, 0, srcg->datasize);
  return 0;
}

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
  if(0 != reset_srcs(fluid->srcux + rk_a, fluid->srcux + rk_b, fluid->srcux + rk_g)){
    return 1;
  }
  if(0 != reset_srcs(fluid->srcuy + rk_a, fluid->srcuy + rk_b, fluid->srcuy + rk_g)){
    return 1;
  }
#if NDIMS == 3
  if(0 != reset_srcs(fluid->srcuz + rk_a, fluid->srcuz + rk_b, fluid->srcuz + rk_g)){
    return 1;
  }
#endif
  if(0 != reset_srcs(fluid->srct  + rk_a, fluid->srct  + rk_b, fluid->srct  + rk_g)){
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

