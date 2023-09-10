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
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  // set boundary values | 13
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    P(      0, j) = P(    1, j); // Neumann
    P(isize+1, j) = P(isize, j); // Neumann
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      P(      0, j, k) = P(    1, j, k); // Neumann
      P(isize+1, j, k) = P(isize, j, k); // Neumann
    }
  }
#endif
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
#if NDIMS == 3
    MPI_DOUBLE,
#endif
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, p)){
    return 1;
  }
#if NDIMS == 3
  if(0 != halo_communicate_in_z(domain, dtypes + 1, p)){
    return 1;
  }
#endif
  assign_boundary_conditions_in_x(domain, p->data);
  return 0;
}

