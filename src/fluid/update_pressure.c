#include "param.h"
#include "runge_kutta.h"
#include "memory.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "array_macros/domain/dxf.h"
#include "array_macros/domain/dxc.h"
#include "array_macros/fluid/p.h"
#include "array_macros/fluid/psi.h"

#define BEGIN \
  for(int k = 1; k <= ksize; k++){ \
    for(int j = 1; j <= jsize; j++){ \
      for(int i = 1; i <= isize; i++){
#define END \
      } \
    } \
  }

static inline int add_explicit(
    const domain_t * domain,
    const fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict psi = fluid->psi.data;
  double * restrict p = fluid->p.data;
  BEGIN
    // explicit contribution
    P(i, j, k) += PSI(i, j, k);
  END
  return 0;
}

static inline int add_implicit_x(
    const domain_t * domain,
    const double prefactor,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict dxf = domain->dxf;
  const double * restrict dxc = domain->dxc;
  const double * restrict psi = fluid->psi.data;
  double * restrict p = fluid->p.data;
  BEGIN
    // x implicit contribution
    const double dpsidx_xm = (- PSI(i-1, j  , k  ) + PSI(i  , j  , k  )) / DXC(i  );
    const double dpsidx_xp = (- PSI(i  , j  , k  ) + PSI(i+1, j  , k  )) / DXC(i+1);
    P(i, j, k) -= prefactor / DXF(i  ) * (
        - dpsidx_xm
        + dpsidx_xp
    );
  END
  return 0;
}

static inline int add_implicit_y(
    const domain_t * domain,
    const double prefactor,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double dy = domain->dy;
  const double * restrict psi = fluid->psi.data;
  double * restrict p = fluid->p.data;
  BEGIN
    // y implicit contribution
    const double dpsidy_ym = (- PSI(i  , j-1, k  ) + PSI(i  , j  , k  )) / dy;
    const double dpsidy_yp = (- PSI(i  , j  , k  ) + PSI(i  , j+1, k  )) / dy;
    P(i, j, k) -= prefactor / dy * (
        - dpsidy_ym
        + dpsidy_yp
    );
  END
  return 0;
}

static inline int add_implicit_z(
    const domain_t * domain,
    const double prefactor,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double dz = domain->dz;
  const double * restrict psi = fluid->psi.data;
  double * restrict p = fluid->p.data;
  BEGIN
    // z implicit contribution
    const double dpsidz_zm = (- PSI(i  , j  , k-1) + PSI(i  , j  , k  )) / dz;
    const double dpsidz_zp = (- PSI(i  , j  , k  ) + PSI(i  , j  , k+1)) / dz;
    P(i, j, k) -= prefactor / dz * (
        - dpsidz_zm
        + dpsidz_zp
    );
  END
  return 0;
}

/**
 * @brief update pressure using scalar potential psi
 * @param[in]     domain : information related to MPI domain decomposition
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : scalar potential (in), pressure (out)
 * @return               : error code
 */
int fluid_update_pressure(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
){
  // explicit contribution, always present
  add_explicit(domain, fluid);
  // gamma dt diffusivity / 2
  const double prefactor =
    0.5 * rkcoefs[rkstep][rk_g] * dt * fluid->m_dif;
  // additional corrections if diffusive terms
  //   in the direction is treated implicitly
  if(param_m_implicit_x){
    add_implicit_x(domain, prefactor, fluid);
  }
  if(param_m_implicit_y){
    add_implicit_y(domain, prefactor, fluid);
  }
  if(param_m_implicit_z){
    add_implicit_z(domain, prefactor, fluid);
  }
  // impose boundary conditions and communicate halo cells
  if(0 != fluid_update_boundaries_p(domain, &fluid->p)){
    return 1;
  }
  return 0;
}

