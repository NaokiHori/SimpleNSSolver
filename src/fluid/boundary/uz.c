#include "param.h"
#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid_solver.h"
#include "array_macros/fluid/uz.h"

static int assign_boundary_conditions_in_x(
    const domain_t * domain,
    double * uz
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  // set boundary values
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      UZ(      0, j, k) = param_uz_xm; // no-slip
      UZ(isize+1, j, k) = param_uz_xp; // no-slip
    }
  }
  return 0;
}

/**
 * @brief update boundary values of z velocity
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] uz     : z velocity
 * @return               : error code
 */
int fluid_update_boundaries_uz(
    const domain_t * domain,
    array_t * uz
){
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, uz)){
    return 1;
  }
  if(0 != halo_communicate_in_z(domain, dtypes + 1, uz)){
    return 1;
  }
  assign_boundary_conditions_in_x(domain, uz->data);
  return 0;
}
