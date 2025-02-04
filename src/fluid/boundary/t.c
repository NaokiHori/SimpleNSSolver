#include "param.h"
#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid_solver.h"
#include "array_macros/fluid/t.h"

static int assign_boundary_conditions_in_x(
    const domain_t * domain,
    double * t
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  // set boundary values
  for(int j = 1; j <= jsize; j++){
    T(      0, j) = param_t_xm;
    T(isize+1, j) = param_t_xp;
  }
  return 0;
}

/**
 * @brief update boundary values of temperature
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] t      : temperature
 * @return               : error code
 */
int fluid_update_boundaries_t(
    const domain_t * domain,
    array_t * t
){
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, t)){
    return 1;
  }
  assign_boundary_conditions_in_x(domain, t->data);
  return 0;
}

