#include "param.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "internal.h"

/**
 * @brief couple external forces
 * @param[in]     domain : information related to MPI domain decomposition
 * @param[in,out] fluid  : buoyancy force (in), x RK source term (in/out)
 * @return               : error code
 */
int fluid_couple_external_force(
    const domain_t * domain,
    fluid_t * fluid
){
  // add buoyancy in x when spcified
  if(param_add_buoyancy){
    if(0 != buoyancy_ux(domain, fluid)){
      return 1;
    }
  }
  return 0;
}

