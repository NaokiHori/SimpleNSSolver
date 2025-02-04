#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "internal.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/psi.h"

#define BEGIN \
  for (int j = 1; j <= jsize; j++) { \
    for (int i = 2; i <= isize; i++) {
#define END \
    } \
  }

// correct x velocity
int fluid_correct_velocity_ux (
    const domain_t * domain,
    const double prefactor,
    fluid_t * fluid
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict psi = fluid->psi.data;
  double * restrict ux = fluid->ux.data;
  BEGIN
    const double hx = HXXF(i  );
    const double psi_xm = PSI(i-1, j  );
    const double psi_xp = PSI(i  , j  );
    double * vel = &UX(i, j);
    *vel -= prefactor / hx * (
        - psi_xm
        + psi_xp
    );
  END
  // update boundary and halo cells
  if (0 != fluid_update_boundaries_ux(domain, &fluid->ux)) {
    return 1;
  }
  return 0;
}

