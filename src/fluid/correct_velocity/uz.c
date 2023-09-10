#if NDIMS == 3
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "internal.h"
#include "array_macros/fluid/uz.h"
#include "array_macros/fluid/psi.h"

/**
 * @brief correct uz using scalar potential psi
 * @param[in]     domain    : information about domain decomposition and size
 * @param[in]     prefactor : pre-factor in front of grad psi
 * @param[in,out] fluid     : scalar potential psi (in), uz (out)
 * @return                  : error code
 */
int fluid_correct_velocity_uz(
    const domain_t * domain,
    const double prefactor,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double dz = domain->dz;
  const double * restrict psi = fluid->psi.data;
  double * restrict uz = fluid->uz.data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // correct z velocity | 6
        const double psi_zm = PSI(i  , j  , k-1);
        const double psi_zp = PSI(i  , j  , k  );
        UZ(i, j, k) -= prefactor / dz * (
            + psi_zp
            - psi_zm
        );
      }
    }
  }
  // update boundary and halo cells | 3
  if(0 != fluid_update_boundaries_uz(domain, &fluid->uz)){
    return 1;
  }
  return 0;
}
#endif
