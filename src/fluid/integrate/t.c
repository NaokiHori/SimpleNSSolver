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
#if NDIMS == 3
#include "array_macros/fluid/uz.h"
#endif
#include "array_macros/fluid/t.h"

// store approximation of laplacian
typedef double laplacian_t[3];

typedef struct {
  bool is_initialised;
  laplacian_t * lapx;
  laplacian_t lapy;
#if NDIMS == 3
  laplacian_t lapz;
#endif
} laplacians_t;

static laplacians_t laplacians = {
  .is_initialised = false,
};

// [1 : isize]
#define LAPX(I) lapx[(I)-1]

static int init_lap(
    const domain_t * domain
){
  // Laplacian w.r.t. temp in x | 14
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
  // Laplacian in y | 6
  {
    const double dy = domain->dy;
    laplacians.lapy[0] = + 1. / dy / dy;
    laplacians.lapy[1] = - 2. / dy / dy;
    laplacians.lapy[2] = + 1. / dy / dy;
  }
#if NDIMS == 3
  // Laplacian in z | 6
  {
    const double dz = domain->dz;
    laplacians.lapz[0] = + 1. / dz / dz;
    laplacians.lapz[1] = - 2. / dz / dz;
    laplacians.lapz[2] = + 1. / dz / dz;
  }
#endif
  laplacians.is_initialised = true;
  return 0;
}

#if NDIMS == 2
#define BEGIN \
  for(int cnt = 0, j = 1; j <= jsize; j++){ \
    for(int i = 1; i <= isize; i++, cnt++){
#define END \
    } \
  }
#else
#define BEGIN \
  for(int cnt = 0, k = 1; k <= ksize; k++){ \
    for(int j = 1; j <= jsize; j++){ \
      for(int i = 1; i <= isize; i++, cnt++){
#define END \
      } \
    } \
  }
#endif

static int advection_x(
    const domain_t * domain,
    const double * restrict t,
    const double * restrict ux,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict dxf = domain->dxf;
#if NDIMS == 2
  BEGIN
    // T is transported by ux | 7
    const double l = + 0.5 / DXF(i  ) * UX(i  , j  );
    const double u = - 0.5 / DXF(i  ) * UX(i+1, j  );
    const double c = - l - u;
    src[cnt] +=
      + l * T(i-1, j  )
      + c * T(i  , j  )
      + u * T(i+1, j  );
  END
#else
  BEGIN
    // T is transported by ux | 7
    const double l = + 0.5 / DXF(i  ) * UX(i  , j  , k  );
    const double u = - 0.5 / DXF(i  ) * UX(i+1, j  , k  );
    const double c = - l - u;
    src[cnt] +=
      + l * T(i-1, j  , k  )
      + c * T(i  , j  , k  )
      + u * T(i+1, j  , k  );
  END
#endif
  return 0;
}

static int advection_y(
    const domain_t * domain,
    const double * restrict t,
    const double * restrict uy,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double dy = domain->dy;
#if NDIMS == 2
  BEGIN
    // T is transported by uy | 7
    const double l = + 0.5 / dy * UY(i  , j  );
    const double u = - 0.5 / dy * UY(i  , j+1);
    const double c = - l - u;
    src[cnt] +=
      + l * T(i  , j-1)
      + c * T(i  , j  )
      + u * T(i  , j+1);
  END
#else
  BEGIN
    // T is transported by uy | 7
    const double l = + 0.5 / dy * UY(i  , j  , k  );
    const double u = - 0.5 / dy * UY(i  , j+1, k  );
    const double c = - l - u;
    src[cnt] +=
      + l * T(i  , j-1, k  )
      + c * T(i  , j  , k  )
      + u * T(i  , j+1, k  );
  END
#endif
  return 0;
}

#if NDIMS == 3
static int advection_z(
    const domain_t * domain,
    const double * restrict t,
    const double * restrict uz,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double dz = domain->dz;
  BEGIN
    // T is transported by uz | 7
    const double l = + 0.5 / dz * UZ(i  , j  , k  );
    const double u = - 0.5 / dz * UZ(i  , j  , k+1);
    const double c = - l - u;
    src[cnt] +=
      + l * T(i  , j  , k-1)
      + c * T(i  , j  , k  )
      + u * T(i  , j  , k+1);
  END
  return 0;
}
#endif

static int diffusion_x(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict t,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const laplacian_t * restrict lapx = laplacians.lapx;
#if NDIMS == 2
  BEGIN
    // T is diffused in x | 5
    src[cnt] += diffusivity * (
        + LAPX(i)[0] * T(i-1, j  )
        + LAPX(i)[1] * T(i  , j  )
        + LAPX(i)[2] * T(i+1, j  )
    );
  END
#else
  BEGIN
    // T is diffused in x | 5
    src[cnt] += diffusivity * (
        + LAPX(i)[0] * T(i-1, j  , k  )
        + LAPX(i)[1] * T(i  , j  , k  )
        + LAPX(i)[2] * T(i+1, j  , k  )
    );
  END
#endif
  return 0;
}

static int diffusion_y(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict t,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const laplacian_t * restrict lapy = &laplacians.lapy;
#if NDIMS == 2
  BEGIN
    // T is diffused in y | 5
    src[cnt] += diffusivity * (
        + (*lapy)[0] * T(i  , j-1)
        + (*lapy)[1] * T(i  , j  )
        + (*lapy)[2] * T(i  , j+1)
    );
  END
#else
  BEGIN
    // T is diffused in y | 5
    src[cnt] += diffusivity * (
        + (*lapy)[0] * T(i  , j-1, k  )
        + (*lapy)[1] * T(i  , j  , k  )
        + (*lapy)[2] * T(i  , j+1, k  )
    );
  END
#endif
  return 0;
}

#if NDIMS == 3
static int diffusion_z(
    const domain_t * domain,
    const double diffusivity,
    const double * restrict t,
    double * restrict src
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t * restrict lapz = &laplacians.lapz;
  BEGIN
    // T is diffused in z | 5
    src[cnt] += diffusivity * (
        + (*lapz)[0] * T(i  , j  , k-1)
        + (*lapz)[1] * T(i  , j  , k  )
        + (*lapz)[2] * T(i  , j  , k+1)
    );
  END
  return 0;
}
#endif

/**
 * @brief comute right-hand-side of Runge-Kutta scheme
 * @param[in]     domain : information related to domain decomposition and size
 * @param[in,out] fluid  : n-step flow field (in), RK source terms (out)
 * @return               : error code
 */
int compute_rhs_t(
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
#if NDIMS == 3
  const double * restrict uz = fluid->uz.data;
#endif
  const double * restrict  t = fluid-> t.data;
  double * restrict srca = fluid->srct[rk_a].data;
  double * restrict srcg = fluid->srct[rk_g].data;
  const double diffusivity = fluid->t_dif;
  // advective contributions, always explicit | 5
  advection_x(domain, t, ux, srca);
  advection_y(domain, t, uy, srca);
#if NDIMS == 3
  advection_z(domain, t, uz, srca);
#endif
  // diffusive contributions, can be explicit or implicit | 5
  diffusion_x(domain, diffusivity, t, param_t_implicit_x ? srcg : srca);
  diffusion_y(domain, diffusivity, t, param_t_implicit_y ? srcg : srca);
#if NDIMS == 3
  diffusion_z(domain, diffusivity, t, param_t_implicit_z ? srcg : srca);
#endif
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

#if NDIMS == 3
static int solve_in_z(
    const double prefactor,
    linear_system_t * linear_system
){
  tdm_info_t * tdm_info = linear_system->tdm_z;
  int size = 0;
  double * restrict tdm_l = NULL;
  double * restrict tdm_c = NULL;
  double * restrict tdm_u = NULL;
  tdm.get_size(tdm_info, &size);
  tdm.get_l(tdm_info, &tdm_l);
  tdm.get_c(tdm_info, &tdm_c);
  tdm.get_u(tdm_info, &tdm_u);
  const laplacian_t * restrict lapz = &laplacians.lapz;
  for(int k = 0; k < size; k++){
    tdm_l[k] =    - prefactor * (*lapz)[0];
    tdm_c[k] = 1. - prefactor * (*lapz)[1];
    tdm_u[k] =    - prefactor * (*lapz)[2];
  }
  tdm.solve(tdm_info, linear_system->z2pncl);
  return 0;
}
#endif

/**
 * @brief update temperature field
 * @param[in]     domain : information about domain decomposition and size
 * @param[in]     rkstep : Runge-Kutta step
 * @param[in]     dt     : time step size
 * @param[in,out] fluid  : Runge-Kutta source terms (in), temperature (out)
 * @return               : error code
 */
int update_t(
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
      param_t_implicit_x,
      param_t_implicit_y,
#if NDIMS == 3
      param_t_implicit_z,
#endif
    };
    const size_t glsizes[NDIMS] = {
      domain->glsizes[0],
      domain->glsizes[1],
#if NDIMS == 3
      domain->glsizes[2],
#endif
    };
    if(0 != linear_system_init(domain->info, implicit, glsizes, &linear_system)){
      return 1;
    }
  }
  // compute increments | 25
  {
    const double coef_a = rkcoefs[rkstep][rk_a];
    const double coef_b = rkcoefs[rkstep][rk_b];
    const double coef_g = rkcoefs[rkstep][rk_g];
    const double * restrict srcta = fluid->srct[rk_a].data;
    const double * restrict srctb = fluid->srct[rk_b].data;
    const double * restrict srctg = fluid->srct[rk_g].data;
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
#if NDIMS == 3
    const int ksize = domain->mysizes[2];
#endif
    double * restrict dtemp = linear_system.x1pncl;
#if NDIMS == 2
    const size_t nitems = isize * jsize;
#else
    const size_t nitems = isize * jsize * ksize;
#endif
    for(size_t n = 0; n < nitems; n++){
      dtemp[n] =
        + coef_a * dt * srcta[n]
        + coef_b * dt * srctb[n]
        + coef_g * dt * srctg[n];
    }
  }
  // gamma dt diffusivity / 2 | 2
  const double prefactor =
    0.5 * rkcoefs[rkstep][rk_g] * dt * fluid->t_dif;
  // solve linear systems in x
  if(param_t_implicit_x){
    solve_in_x(
        prefactor,
        &linear_system
    );
  }
  // solve linear systems in y
  if(param_t_implicit_y){
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
#if NDIMS == 3
  // solve linear systems in z
  if(param_t_implicit_z){
    sdecomp.transpose.execute(
        linear_system.transposer_x1_to_z2,
        linear_system.x1pncl,
        linear_system.z2pncl
    );
    solve_in_z(
        prefactor,
        &linear_system
    );
    sdecomp.transpose.execute(
        linear_system.transposer_z2_to_x1,
        linear_system.z2pncl,
        linear_system.x1pncl
    );
  }
#endif
  // the field is actually updated here | 21
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
#if NDIMS == 3
    const int ksize = domain->mysizes[2];
#endif
    const double * restrict dtemp = linear_system.x1pncl;
    double * restrict t = fluid->t.data;
#if NDIMS == 2
    BEGIN
      T(i, j) += dtemp[cnt];
    END
#else
    BEGIN
      T(i, j, k) += dtemp[cnt];
    END
#endif
    if(0 != fluid_update_boundaries_t(domain, &fluid->t)){
      return 1;
    }
  }
  return 0;
}

