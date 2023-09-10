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

#if NDIMS == 2
#define BEGIN \
  for(int j = 1; j <= jsize; j++){ \
    for(int i = 1; i <= isize; i++){
#define END \
    } \
  }
#else
#define BEGIN \
  for(int k = 1; k <= ksize; k++){ \
    for(int j = 1; j <= jsize; j++){ \
      for(int i = 1; i <= isize; i++){
#define END \
      } \
    } \
  }
#endif

static inline int add_explicit(
    const domain_t * domain,
    const fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict psi = fluid->psi.data;
  double * restrict p = fluid->p.data;
#if NDIMS == 2
  BEGIN
    // explicit contribution | 1
    P(i, j) += PSI(i, j);
  END
#else
  BEGIN
    // explicit contribution | 1
    P(i, j, k) += PSI(i, j, k);
  END
#endif
  return 0;
}

static inline int add_implicit_x(
    const domain_t * domain,
    const double prefactor,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict dxf = domain->dxf;
  const double * restrict dxc = domain->dxc;
  const double * restrict psi = fluid->psi.data;
  double * restrict p = fluid->p.data;
#if NDIMS == 2
  BEGIN
    // x implicit contribution | 6
    const double dpsidx_xm = (- PSI(i-1, j  ) + PSI(i  , j  )) / DXC(i  );
    const double dpsidx_xp = (- PSI(i  , j  ) + PSI(i+1, j  )) / DXC(i+1);
    P(i, j) -= prefactor / DXF(i  ) * (
        - dpsidx_xm
        + dpsidx_xp
    );
  END
#else
  BEGIN
    // x implicit contribution | 6
    const double dpsidx_xm = (- PSI(i-1, j  , k  ) + PSI(i  , j  , k  )) / DXC(i  );
    const double dpsidx_xp = (- PSI(i  , j  , k  ) + PSI(i+1, j  , k  )) / DXC(i+1);
    P(i, j, k) -= prefactor / DXF(i  ) * (
        - dpsidx_xm
        + dpsidx_xp
    );
  END
#endif
  return 0;
}

static inline int add_implicit_y(
    const domain_t * domain,
    const double prefactor,
    fluid_t * fluid
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double dy = domain->dy;
  const double * restrict psi = fluid->psi.data;
  double * restrict p = fluid->p.data;
#if NDIMS == 2
  BEGIN
    // y implicit contribution | 6
    const double dpsidy_ym = (- PSI(i  , j-1) + PSI(i  , j  )) / dy;
    const double dpsidy_yp = (- PSI(i  , j  ) + PSI(i  , j+1)) / dy;
    P(i, j) -= prefactor / dy * (
        - dpsidy_ym
        + dpsidy_yp
    );
  END
#else
  BEGIN
    // y implicit contribution | 6
    const double dpsidy_ym = (- PSI(i  , j-1, k  ) + PSI(i  , j  , k  )) / dy;
    const double dpsidy_yp = (- PSI(i  , j  , k  ) + PSI(i  , j+1, k  )) / dy;
    P(i, j, k) -= prefactor / dy * (
        - dpsidy_ym
        + dpsidy_yp
    );
  END
#endif
  return 0;
}

#if NDIMS == 3
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
    // z implicit contribution | 6
    const double dpsidz_zm = (- PSI(i  , j  , k-1) + PSI(i  , j  , k  )) / dz;
    const double dpsidz_zp = (- PSI(i  , j  , k  ) + PSI(i  , j  , k+1)) / dz;
    P(i, j, k) -= prefactor / dz * (
        - dpsidz_zm
        + dpsidz_zp
    );
  END
  return 0;
}
#endif

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
  // gamma dt diffusivity / 2 | 2
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
#if NDIMS == 3
  if(param_m_implicit_z){
    add_implicit_z(domain, prefactor, fluid);
  }
#endif
  // impose boundary conditions and communicate halo cells | 3
  if(0 != fluid_update_boundaries_p(domain, &fluid->p)){
    return 1;
  }
  return 0;
}

