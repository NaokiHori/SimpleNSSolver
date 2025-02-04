#include "param.h"
#include "runge_kutta.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/p.h"
#include "array_macros/fluid/psi.h"

#define BEGIN \
  for (int j = 1; j <= jsize; j++) { \
    for (int i = 1; i <= isize; i++) {
#define END \
    } \
  }

// explicit contribution
static int explicit_contribution (
    const domain_t * domain,
    const fluid_t * fluid
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict psi = fluid->psi.data;
  double * restrict p = fluid->p.data;
  BEGIN
    P(i, j) += PSI(i, j);
  END
  return 0;
}

// x implicit contribution
static int implicit_x_contribution (
    const domain_t * domain,
    const double prefactor,
    fluid_t * fluid
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict psi = fluid->psi.data;
  double * restrict p = fluid->p.data;
  BEGIN
    const double hx_xm = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i  );
    const double jd_x0 = JDXC(i  );
    const double jd_xp = JDXF(i+1);
    const double l = 1. / jd_x0 * jd_xm / hx_xm / hx_xm;
    const double u = 1. / jd_x0 * jd_xp / hx_xp / hx_xp;
    const double c = - l - u;
    const double psi_xm = PSI(i-1, j  );
    const double psi_x0 = PSI(i  , j  );
    const double psi_xp = PSI(i+1, j  );
    double * pre = &P(i, j);
    *pre -= prefactor * (
        + l * psi_xm
        + c * psi_x0
        + u * psi_xp
    );
  END
  return 0;
}

// y implicit contribution
static int implicit_y_contribution (
    const domain_t * domain,
    const double prefactor,
    fluid_t * fluid
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double hy = domain->hy;
  const double * restrict psi = fluid->psi.data;
  double * restrict p = fluid->p.data;
  BEGIN
    const double l = 1. / hy / hy;
    const double u = 1. / hy / hy;
    const double c = - l - u;
    const double psi_ym = PSI(i  , j-1);
    const double psi_y0 = PSI(i  , j  );
    const double psi_yp = PSI(i  , j+1);
    double * pre = &P(i, j);
    *pre -= prefactor * (
        + l * psi_ym
        + c * psi_y0
        + u * psi_yp
    );
  END
  return 0;
}

// update pressure field using scalar potential
int fluid_update_pressure (
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
) {
  // explicit contribution, always present
  explicit_contribution(domain, fluid);
  const double prefactor =
    0.5 * rkcoefs[rkstep].gamma * dt * fluid_compute_momentum_diffusivity(fluid);
  // additional terms if diffusive terms in the direction is treated implicitly
  if (param_m_implicit_x) {
    implicit_x_contribution(domain, prefactor, fluid);
  }
  if (param_m_implicit_y) {
    implicit_y_contribution(domain, prefactor, fluid);
  }
  // impose boundary conditions and communicate halo cells
  if (0 != fluid_update_boundaries_p(domain, &fluid->p)) {
    return 1;
  }
  return 0;
}

