#include <string.h>
#include "runge_kutta.h"
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
 * @brief compute right-hand-side terms of the Runge-Kutta scheme
 * @param[in]     domain : information related to MPI domain decomposition
 * @param[in,out] fluid  : flow field (in), RK source terms (out)
 * @return               : error code
 */
int fluid_compute_rhs(
    const domain_t * domain,
    fluid_t * fluid
){
  // reset buffers
  // copy previous k-step source term and reset
  if(0 != reset_srcs(fluid->srcux + rk_a, fluid->srcux + rk_b, fluid->srcux + rk_g)){
    return 1;
  }
  if(0 != reset_srcs(fluid->srcuy + rk_a, fluid->srcuy + rk_b, fluid->srcuy + rk_g)){
    return 1;
  }
  if(0 != reset_srcs(fluid->srct  + rk_a, fluid->srct  + rk_b, fluid->srct  + rk_g)){
    return 1;
  }
  // update buffers
  if(0 != compute_rhs_ux(domain, fluid)){
    return 1;
  }
  if(0 != compute_rhs_uy(domain, fluid)){
    return 1;
  }
  if(0 != compute_rhs_t (domain, fluid)){
    return 1;
  }
  return 0;
}

