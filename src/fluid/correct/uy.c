#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "internal.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/psi.h"

#define BEGIN \
  for (int j = 1; j <= jsize; j++) { \
    for (int i = 1; i <= isize; i++) {
#define END \
    } \
  }

// correct y velocity
int fluid_correct_velocity_uy (
    const domain_t * domain,
    const double prefactor,
    fluid_t * fluid
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double hy = domain->hy;
  const double * restrict psi = fluid->psi.data;
  double * restrict uy = fluid->uy.data;
  BEGIN
    const double psi_ym = PSI(i  , j-1);
    const double psi_yp = PSI(i  , j  );
    double * vel = &UY(i, j);
    *vel -= prefactor / hy * (
        - psi_ym
        + psi_yp
    );
  END
  // update boundary and halo cells
  if (0 != fluid_update_boundaries_uy(domain, &fluid->uy)) {
    return 1;
  }
  return 0;
}

