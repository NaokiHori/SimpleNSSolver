#include "array.h"
#include "domain.h"
#include "halo.h"
#include "fluid_solver.h"
#include "array_macros/fluid/ux.h"

static int assign_boundary_conditions_in_x(
    const domain_t * domain,
    double * ux
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  // set boundary values | 13
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    UX(      1, j) = 0.; // impermeable
    UX(isize+1, j) = 0.; // impermeable
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      UX(      1, j, k) = 0.; // impermeable
      UX(isize+1, j, k) = 0.; // impermeable
    }
  }
#endif
  return 0;
}

/**
 * @brief update boundary values of x velocity
 * @param[in]     domain : information about domain decomposition and size
 * @param[in,out] ux     : x velocity
 * @return               : error code
 */
int fluid_update_boundaries_ux(
    const domain_t * domain,
    array_t * ux
){
  static MPI_Datatype dtypes[NDIMS - 1] = {
    MPI_DOUBLE,
#if NDIMS == 3
    MPI_DOUBLE,
#endif
  };
  if(0 != halo_communicate_in_y(domain, dtypes + 0, ux)){
    return 1;
  }
#if NDIMS == 3
  if(0 != halo_communicate_in_z(domain, dtypes + 1, ux)){
    return 1;
  }
#endif
  assign_boundary_conditions_in_x(domain, ux->data);
  return 0;
}

