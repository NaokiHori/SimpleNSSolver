#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "internal.h"
#include "array_macros/domain/dxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/psi.h"

/**
 * @brief correct ux using scalar potential psi
 * @param[in]     domain    : information about domain decomposition and size
 * @param[in]     prefactor : pre-factor in front of grad psi
 * @param[in,out] fluid     : scalar potential psi (in), ux (out)
 * @return                  : error code
 */
int fluid_correct_velocity_ux(
    const domain_t * domain,
    const double prefactor,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict dxc = domain->dxc;
  const double * restrict psi = fluid->psi.data;
  double * restrict ux = fluid->ux.data;
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      // correct x velocity | 7
      const double dx = DXC(i  );
      const double psi_xm = PSI(i-1, j  );
      const double psi_xp = PSI(i  , j  );
      UX(i, j) -= prefactor / dx * (
          + psi_xp
          - psi_xm
      );
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        // correct x velocity | 7
        const double dx = DXC(i  );
        const double psi_xm = PSI(i-1, j  , k  );
        const double psi_xp = PSI(i  , j  , k  );
        UX(i, j, k) -= prefactor / dx * (
            + psi_xp
            - psi_xm
        );
      }
    }
  }
#endif
  // update boundary and halo cells | 3
  if(0 != fluid_update_boundaries_ux(domain, &fluid->ux)){
    return 1;
  }
  return 0;
}

