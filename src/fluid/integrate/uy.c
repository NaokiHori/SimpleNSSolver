#include "param.h"
#include "memory.h"
#include "runge_kutta.h"
#include "linear_system.h"
#include "tdm.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "internal.h"
#include "array_macros/domain/dxf.h"
#include "array_macros/domain/dxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/p.h"

// store approximation of laplacian
typedef double laplacian_t[3];

typedef struct {
  bool is_initialised;
  laplacian_t * lapx;
  laplacian_t lapy;
} laplacians_t;

static laplacians_t laplacians = {
  .is_initialised = false,
};

// [1 : isize]
#define LAPX(I) lapx[(I)-1]

static int init_lap(
    const domain_t * domain
){
  // Laplacian w.r.t. uy in x
  {
    const size_t isize = domain->glsizes[0];
    const double * dxf = domain->dxf;
    const double * dxc = domain->dxc;
    laplacians.lapx = memory_calloc(isize, sizeof(laplacian_t));
    for(size_t i = 1; i <= isize; i++){
      const double l = 1. / DXC(i  ) / DXF(i  );
      const double u = 1. / DXC(i+1) / DXF(i  );
      const double c = - l - u;
      laplacians.LAPX(i)[0] = l;
      laplacians.LAPX(i)[1] = c;
      laplacians.LAPX(i)[2] = u;
    }
  }
  // Laplacian in y
  {
    const double dy = domain->dy;
    laplacians.lapy[0] = + 1. / dy / dy;
    laplacians.lapy[1] = - 2. / dy / dy;
    laplacians.lapy[2] = + 1. / dy / dy;
  }
  laplacians.is_initialised = true;
  return 0;
}

#define BEGIN \
  for(int cnt = 0, j = 1; j <= jsize; j++){ \
    for(int i = 1; i <= isize; i++, cnt++){
#define END \
    } \
  }

static int advection_x(
    const domain_t * domain,
    const double * restrict uy,
    const double * restrict ux,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict dxf = domain->dxf;
  BEGIN
    // uy is transported by ux
    const double ux_l = + 0.5 * UX(i  , j-1) + 0.5 * UX(i  , j  );
    const double ux_u = + 0.5 * UX(i+1, j-1) + 0.5 * UX(i+1, j  );
    const double l = + 0.5 / DXF(i  ) * ux_l;
    const double u = - 0.5 / DXF(i  ) * ux_u;
    const double c = - l - u;
    src[cnt] +=
      + l * UY(i-1, j  )
      + c * UY(i  , j  )
      + u * UY(i+1, j  );
  END
  return 0;
}

static int advection_y(
    const domain_t * domain,
    const double * restrict uy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double dy = domain->dy;
  BEGIN
    // uy is transported by uy
    const double uy_l = + 0.5 * UY(i  , j-1) + 0.5 * UY(i  , j  );
    const double uy_u = + 0.5 * UY(i  , j  ) + 0.5 * UY(i  , j+1);
    const double l = + 0.5 / dy * uy_l;
    const double u = - 0.5 / dy * uy_u;
    const double c = - l - u;
    src[cnt] +=
      + l * UY(i  , j-1)
      + c * UY(i  , j  )
      + u * UY(i  , j+1);
  END
  return 0;
}

static int diffusion_x(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict uy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const laplacian_t * restrict lapx = laplacians.lapx;
  BEGIN
    // uy is diffused in x
    src[cnt] += diffusivity * (
        + LAPX(i)[0] * UY(i-1, j  )
        + LAPX(i)[1] * UY(i  , j  )
        + LAPX(i)[2] * UY(i+1, j  )
    );
  END
  return 0;
}

static int diffusion_y(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict uy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const laplacian_t * restrict lapy = &laplacians.lapy;
  BEGIN
    // uy is diffused in y
    src[cnt] += diffusivity * (
        + (*lapy)[0] * UY(i  , j-1)
        + (*lapy)[1] * UY(i  , j  )
        + (*lapy)[2] * UY(i  , j+1)
    );
  END
  return 0;
}

static int pressure(
    const domain_t * domain,
    const double * restrict p,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double dy = domain->dy;
  BEGIN
    src[cnt] -= 1. / dy * (
        - P(i  , j-1)
        + P(i  , j  )
    );
  END
  return 0;
}

/**
 * @brief comute right-hand-side of Runge-Kutta scheme of uy
 * @param[in]     domain : information related to MPI domain decomposition
 * @param[in,out] fluid  : n-step flow field (in), RK source terms (inout)
 * @return               : error code
 */
int compute_rhs_uy(
    const domain_t * domain,
    fluid_t * fluid
){
  if(!laplacians.is_initialised){
    if(0 != init_lap(domain)){
      return 1;
    }
  }
  const double * restrict ux = fluid->ux.data;
  const double * restrict uy = fluid->uy.data;
  const double * restrict  p = fluid-> p.data;
  double * restrict srca = fluid->srcuy[rk_a].data;
  double * restrict srcg = fluid->srcuy[rk_g].data;
  const double diffusivity = fluid->m_dif;
  // advective contributions, always explicit
  advection_x(domain, uy, ux, srca);
  advection_y(domain, uy,     srca);
  // diffusive contributions, can be explicit or implicit
  diffusion_x(domain, diffusivity, uy, param_m_implicit_x ? srcg : srca);
  diffusion_y(domain, diffusivity, uy, param_m_implicit_y ? srcg : srca);
  // pressure-gradient contribution, always implicit
  pressure(domain, p, srcg);
  return 0;
}

static int solve_in_x(
    const double prefactor,
    linear_system_t * linear_system
){
  tdm_info_t * tdm_info = linear_system->tdm_x;
  int size = 0;
  double * restrict tdm_l = NULL;
  double * restrict tdm_c = NULL;
  double * restrict tdm_u = NULL;
  tdm.get_size(tdm_info, &size);
  tdm.get_l(tdm_info, &tdm_l);
  tdm.get_c(tdm_info, &tdm_c);
  tdm.get_u(tdm_info, &tdm_u);
  const laplacian_t * restrict lapx = laplacians.lapx;
  for(int i = 0; i < size; i++){
    tdm_l[i] =    - prefactor * lapx[i][0];
    tdm_c[i] = 1. - prefactor * lapx[i][1];
    tdm_u[i] =    - prefactor * lapx[i][2];
  }
  tdm.solve(tdm_info, linear_system->x1pncl);
  return 0;
}

static int solve_in_y(
    const double prefactor,
    linear_system_t * linear_system
){
  tdm_info_t * tdm_info = linear_system->tdm_y;
  int size = 0;
  double * restrict tdm_l = NULL;
  double * restrict tdm_c = NULL;
  double * restrict tdm_u = NULL;
  tdm.get_size(tdm_info, &size);
  tdm.get_l(tdm_info, &tdm_l);
  tdm.get_c(tdm_info, &tdm_c);
  tdm.get_u(tdm_info, &tdm_u);
  const laplacian_t * restrict lapy = &laplacians.lapy;
  for(int j = 0; j < size; j++){
    tdm_l[j] =    - prefactor * (*lapy)[0];
    tdm_c[j] = 1. - prefactor * (*lapy)[1];
    tdm_u[j] =    - prefactor * (*lapy)[2];
  }
  tdm.solve(tdm_info, linear_system->y1pncl);
  return 0;
}

/**
 * @brief update uy
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : Runge-Kutta source terms (in), velocity (out)
 * @return               : error code
 */
int update_uy(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
){
  static linear_system_t linear_system = {
    .is_initialised = false,
  };
  if(!linear_system.is_initialised){
    // if not initialised yet, prepare linear solver
    //   for implicit diffusive term treatment
    const bool implicit[NDIMS] = {
      param_m_implicit_x,
      param_m_implicit_y,
    };
    const size_t glsizes[NDIMS] = {
      domain->glsizes[0],
      domain->glsizes[1],
    };
    if(0 != linear_system_init(domain->info, implicit, glsizes, &linear_system)){
      return 1;
    }
  }
  // compute increments
  {
    const double coef_a = rkcoefs[rkstep][rk_a];
    const double coef_b = rkcoefs[rkstep][rk_b];
    const double coef_g = rkcoefs[rkstep][rk_g];
    const double * restrict srcuya = fluid->srcuy[rk_a].data;
    const double * restrict srcuyb = fluid->srcuy[rk_b].data;
    const double * restrict srcuyg = fluid->srcuy[rk_g].data;
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    double * restrict duy = linear_system.x1pncl;
    const size_t nitems = isize * jsize;
    for(size_t n = 0; n < nitems; n++){
      duy[n] =
        + coef_a * dt * srcuya[n]
        + coef_b * dt * srcuyb[n]
        + coef_g * dt * srcuyg[n];
    }
  }
  // gamma dt diffusivity / 2
  const double prefactor =
    0.5 * rkcoefs[rkstep][rk_g] * dt * fluid->m_dif;
  // solve linear systems in x
  if(param_m_implicit_x){
    solve_in_x(
        prefactor,
        &linear_system
    );
  }
  // solve linear systems in y
  if(param_m_implicit_y){
    sdecomp.transpose.execute(
        linear_system.transposer_x1_to_y1,
        linear_system.x1pncl,
        linear_system.y1pncl
    );
    solve_in_y(
        prefactor,
        &linear_system
    );
    sdecomp.transpose.execute(
        linear_system.transposer_y1_to_x1,
        linear_system.y1pncl,
        linear_system.x1pncl
    );
  }
  // the field is actually updated here
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const double * restrict duy = linear_system.x1pncl;
    double * restrict uy = fluid->uy.data;
    BEGIN
      UY(i, j) += duy[cnt];
    END
    if(0 != fluid_update_boundaries_uy(domain, &fluid->uy)){
      return 1;
    }
  }
  return 0;
}

