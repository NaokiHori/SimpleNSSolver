#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid_solver.h"
#include "array_macros/fluid/p.h"

static int assign_boundary_conditions_in_x(
    const domain_t * domain,
    double * p
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  // set boundary values
  for(int j = 1; j <= jsize; j++){
    P(      0, j) = P(    1, j); // Neumann
    P(isize+1, j) = P(isize, j); // Neumann
  }
  return 0;
}

/**
 * @brief update boundary values of the pressure
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] p      : pressure
 * @return               : error code
 */
int fluid_update_boundaries_p(
    const domain_t * domain,
    array_t * p
){
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, p)){
    return 1;
  }
  assign_boundary_conditions_in_x(domain, p->data);
  return 0;
}

