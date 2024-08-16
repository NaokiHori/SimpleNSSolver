#include "param.h"
#include "memory.h"
#include "runge_kutta.h"
#include "linear_system.h"
#include "tdm.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "internal.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/hxxc.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#if NDIMS == 3
#include "array_macros/fluid/uz.h"
#endif
#include "array_macros/fluid/t.h"

static laplacians_t laplacians = {
  .is_initialised = false,
};

// map [1 : isize] in code to [0 : isize - 1] in memory
#define LAPX(I) lapx[(I) - 1]

static int init_laplacians (
    const domain_t * domain
) {
  // scalar laplacian in x | 15
  {
    const int isize = domain->glsizes[0];
    const double * hxxf = domain->hxxf;
    const double * jdxf = domain->jdxf;
    const double * jdxc = domain->jdxc;
    laplacians.lapx = memory_calloc(isize, sizeof(laplacian_t));
    for (int i = 1; i <= isize; i++) {
      const double l = 1. / JDXC(i  ) * JDXF(i  ) / HXXF(i  ) / HXXF(i  );
      const double u = 1. / JDXC(i  ) * JDXF(i+1) / HXXF(i+1) / HXXF(i+1);
      const double c = - l - u;
      laplacians.LAPX(i).l = l;
      laplacians.LAPX(i).c = c;
      laplacians.LAPX(i).u = u;
    }
  }
  // scalar laplacian in y | 9
  {
    const double hy = domain->hy;
    const double l = 1. / hy / hy;
    const double u = 1. / hy / hy;
    const double c = - l - u;
    laplacians.lapy.l = l;
    laplacians.lapy.c = c;
    laplacians.lapy.u = u;
  }
#if NDIMS == 3
  // scalar laplacian in z | 9
  {
    const double hz = domain->hz;
    const double l = 1. / hz / hz;
    const double u = 1. / hz / hz;
    const double c = - l - u;
    laplacians.lapz.l = l;
    laplacians.lapz.c = c;
    laplacians.lapz.u = u;
  }
#endif
  laplacians.is_initialised = true;
  return 0;
}

#if NDIMS == 2
#define BEGIN \
  for (int cnt = 0, j = 1; j <= jsize; j++) { \
    for (int i = 1; i <= isize; i++, cnt++) {
#define END \
    } \
  }
#else
#define BEGIN \
  for (int cnt = 0, k = 1; k <= ksize; k++) { \
    for (int j = 1; j <= jsize; j++) { \
      for (int i = 1; i <= isize; i++, cnt++) {
#define END \
      } \
    } \
  }
#endif

// advected in x | 44
static int advection_x (
    const domain_t * domain,
    const double * restrict t,
    const double * restrict ux,
    double * restrict src
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  BEGIN
    const double hx_xm = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i  );
    const double jd_x0 = JDXC(i  );
    const double jd_xp = JDXF(i+1);
#if NDIMS == 2
    const double ux_xm = jd_xm / hx_xm * UX(i  , j  );
    const double ux_xp = jd_xp / hx_xp * UX(i+1, j  );
#else
    const double ux_xm = jd_xm / hx_xm * UX(i  , j  , k  );
    const double ux_xp = jd_xp / hx_xp * UX(i+1, j  , k  );
#endif
    const double l = - 0.5 * ux_xm;
    const double u = + 0.5 * ux_xp;
    const double c = - l - u;
    src[cnt] -= 1. / jd_x0 * (
#if NDIMS == 2
        + l * T(i-1, j  )
        + c * T(i  , j  )
        + u * T(i+1, j  )
#else
        + l * T(i-1, j  , k  )
        + c * T(i  , j  , k  )
        + u * T(i+1, j  , k  )
#endif
    );
  END
  return 0;
}

// advected in y | 39
static int advection_y (
    const domain_t * domain,
    const double * restrict t,
    const double * restrict uy,
    double * restrict src
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double hy = domain->hy;
  const double * restrict jdxc = domain->jdxc;
  BEGIN
    const double jd = JDXC(i  );
#if NDIMS == 2
    const double uy_ym = jd / hy * UY(i  , j  );
    const double uy_yp = jd / hy * UY(i  , j+1);
#else
    const double uy_ym = jd / hy * UY(i  , j  , k  );
    const double uy_yp = jd / hy * UY(i  , j+1, k  );
#endif
    const double l = - 0.5 * uy_ym;
    const double u = + 0.5 * uy_yp;
    const double c = - l - u;
    src[cnt] -= 1. / jd * (
#if NDIMS == 2
        + l * T(i  , j-1)
        + c * T(i  , j  )
        + u * T(i  , j+1)
#else
        + l * T(i  , j-1, k  )
        + c * T(i  , j  , k  )
        + u * T(i  , j+1, k  )
#endif
    );
  END
  return 0;
}

#if NDIMS == 3
// advected in z | 26
static int advection_z (
    const domain_t * domain,
    const double * restrict t,
    const double * restrict uz,
    double * restrict src
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  const double * restrict jdxc = domain->jdxc;
  BEGIN
    const double jd = JDXC(i  );
    const double uz_zm = jd / hz * UZ(i  , j  , k  );
    const double uz_zp = jd / hz * UZ(i  , j  , k+1);
    const double l = - 0.5 * uz_zm;
    const double u = + 0.5 * uz_zp;
    const double c = - l - u;
    src[cnt] -= 1. / jd * (
        + l * T(i  , j  , k-1)
        + c * T(i  , j  , k  )
        + u * T(i  , j  , k+1)
    );
  END
  return 0;
}
#endif

// diffused in x | 27
static int diffusion_x (
    const domain_t * domain,
    const double diffusivity,
    const double * restrict t,
    double * restrict src
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const laplacian_t * lapx = laplacians.lapx;
  BEGIN
    src[cnt] += diffusivity * (
#if NDIMS == 2
        + LAPX(i).l * T(i-1, j  )
        + LAPX(i).c * T(i  , j  )
        + LAPX(i).u * T(i+1, j  )
#else
        + LAPX(i).l * T(i-1, j  , k  )
        + LAPX(i).c * T(i  , j  , k  )
        + LAPX(i).u * T(i+1, j  , k  )
#endif
    );
  END
  return 0;
}

// diffused in y | 27
static int diffusion_y (
    const domain_t * domain,
    const double diffusivity,
    const double * restrict t,
    double * restrict src
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const laplacian_t * lapy = &laplacians.lapy;
  BEGIN
    src[cnt] += diffusivity * (
#if NDIMS == 2
        + (*lapy).l * T(i  , j-1)
        + (*lapy).c * T(i  , j  )
        + (*lapy).u * T(i  , j+1)
#else
        + (*lapy).l * T(i  , j-1, k  )
        + (*lapy).c * T(i  , j  , k  )
        + (*lapy).u * T(i  , j+1, k  )
#endif
    );
  END
  return 0;
}

#if NDIMS == 3
// diffused in z | 19
static int diffusion_z (
    const domain_t * domain,
    const double diffusivity,
    const double * restrict t,
    double * restrict src
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t * lapz = &laplacians.lapz;
  BEGIN
    src[cnt] += diffusivity * (
        + (*lapz).l * T(i  , j  , k-1)
        + (*lapz).c * T(i  , j  , k  )
        + (*lapz).u * T(i  , j  , k+1)
    );
  END
  return 0;
}
#endif

// compute right-hand-side terms, which are added to buffers | 34
int compute_rhs_t (
    const domain_t * domain,
    fluid_t * fluid
) {
  if (!laplacians.is_initialised) {
    if (0 != init_laplacians(domain)) {
      return 1;
    }
  }
  const double * restrict ux = fluid->ux.data;
  const double * restrict uy = fluid->uy.data;
#if NDIMS == 3
  const double * restrict uz = fluid->uz.data;
#endif
  const double * restrict  t = fluid-> t.data;
  // buffer for explicit terms
  double * restrict srca = fluid->srct.alpha.data;
  // buffer for implicit terms
  double * restrict srcg = fluid->srct.gamma.data;
  const double diffusivity = fluid_compute_temperature_diffusivity(fluid);
  // advective contributions, always explicit
  advection_x(domain, t, ux, srca);
  advection_y(domain, t, uy, srca);
#if NDIMS == 3
  advection_z(domain, t, uz, srca);
#endif
  // diffusive contributions, can be explicit or implicit
  diffusion_x(domain, diffusivity, t, param_t_implicit_x ? srcg : srca);
  diffusion_y(domain, diffusivity, t, param_t_implicit_y ? srcg : srca);
#if NDIMS == 3
  diffusion_z(domain, diffusivity, t, param_t_implicit_z ? srcg : srca);
#endif
  return 0;
}

// update temperature field
int update_t (
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
) {
  static linear_system_t linear_system = {
    .is_initialised = false,
  };
  if (!linear_system.is_initialised) {
    // if not initialised yet, prepare linear solver
    //   for implicit diffusive treatment
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
    if (0 != linear_system_init(domain->info, implicit, glsizes, &linear_system)) {
      return 1;
    }
  }
  // compute increments | 28
  {
    // Runge-Kutta coefficients, alpha, beta, gamma
    const double coef_a = rkcoefs[rkstep].alpha;
    const double coef_b = rkcoefs[rkstep].beta ;
    const double coef_g = rkcoefs[rkstep].gamma;
    // Runge-Kutta buffers, alpha, beta, gamma
    const double * restrict srcta = fluid->srct.alpha.data;
    const double * restrict srctb = fluid->srct.beta .data;
    const double * restrict srctg = fluid->srct.gamma.data;
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
    // compute T(new) - T(old)
    for (size_t n = 0; n < nitems; n++) {
      dtemp[n] =
        + coef_a * dt * srcta[n]
        + coef_b * dt * srctb[n]
        + coef_g * dt * srctg[n];
    }
  }
  // solve linear systems if necessary
  {
    const double prefactor =
      0.5 * rkcoefs[rkstep].gamma * dt * fluid_compute_temperature_diffusivity(fluid);
    // solve linear systems in x | 20
    if (param_t_implicit_x) {
      // prepare tri-diagonal coefficients
      tdm_info_t * tdm_info = linear_system.tdm_x;
      int size = 0;
      double * restrict tdm_l = NULL;
      double * restrict tdm_c = NULL;
      double * restrict tdm_u = NULL;
      tdm.get_size(tdm_info, &size);
      tdm.get_l(tdm_info, &tdm_l);
      tdm.get_c(tdm_info, &tdm_c);
      tdm.get_u(tdm_info, &tdm_u);
      const laplacian_t * lapx = laplacians.lapx;
      for (int i = 0; i < size; i++) {
        tdm_l[i] =    - prefactor * lapx[i].l;
        tdm_c[i] = 1. - prefactor * lapx[i].c;
        tdm_u[i] =    - prefactor * lapx[i].u;
      }
      // solve all
      tdm.solve(tdm_info, linear_system.x1pncl);
    }
    // solve linear systems in y | 32
    if (param_t_implicit_y) {
      // make all y data accessible
      sdecomp.transpose.execute(
          linear_system.transposer_x1_to_y1,
          linear_system.x1pncl,
          linear_system.y1pncl
      );
      // prepare tri-diagonal coefficients
      tdm_info_t * tdm_info = linear_system.tdm_y;
      int size = 0;
      double * restrict tdm_l = NULL;
      double * restrict tdm_c = NULL;
      double * restrict tdm_u = NULL;
      tdm.get_size(tdm_info, &size);
      tdm.get_l(tdm_info, &tdm_l);
      tdm.get_c(tdm_info, &tdm_c);
      tdm.get_u(tdm_info, &tdm_u);
      const laplacian_t * lapy = &laplacians.lapy;
      for (int j = 0; j < size; j++) {
        tdm_l[j] =    - prefactor * (*lapy).l;
        tdm_c[j] = 1. - prefactor * (*lapy).c;
        tdm_u[j] =    - prefactor * (*lapy).u;
      }
      // solve all
      tdm.solve(tdm_info, linear_system.y1pncl);
      // recover original memory alignment
      sdecomp.transpose.execute(
          linear_system.transposer_y1_to_x1,
          linear_system.y1pncl,
          linear_system.x1pncl
      );
    }
#if NDIMS == 3
    // solve linear systems in z | 32
    if (param_t_implicit_z) {
      // make all z data accessible
      sdecomp.transpose.execute(
          linear_system.transposer_x1_to_z2,
          linear_system.x1pncl,
          linear_system.z2pncl
      );
      // prepare tri-diagonal coefficients
      tdm_info_t * tdm_info = linear_system.tdm_z;
      int size = 0;
      double * restrict tdm_l = NULL;
      double * restrict tdm_c = NULL;
      double * restrict tdm_u = NULL;
      tdm.get_size(tdm_info, &size);
      tdm.get_l(tdm_info, &tdm_l);
      tdm.get_c(tdm_info, &tdm_c);
      tdm.get_u(tdm_info, &tdm_u);
      const laplacian_t * lapz = &laplacians.lapz;
      for (int k = 0; k < size; k++) {
        tdm_l[k] =    - prefactor * (*lapz).l;
        tdm_c[k] = 1. - prefactor * (*lapz).c;
        tdm_u[k] =    - prefactor * (*lapz).u;
      }
      // solve all
      tdm.solve(tdm_info, linear_system.z2pncl);
      // recover original memory alignment
      sdecomp.transpose.execute(
          linear_system.transposer_z2_to_x1,
          linear_system.z2pncl,
          linear_system.x1pncl
      );
    }
#endif
  }
  // update temperature field | 19
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
#if NDIMS == 3
    const int ksize = domain->mysizes[2];
#endif
    const double * restrict dtemp = linear_system.x1pncl;
    double * restrict t = fluid->t.data;
    BEGIN
#if NDIMS == 2
      T(i, j) += dtemp[cnt];
#else
      T(i, j, k) += dtemp[cnt];
#endif
    END
    if (0 != fluid_update_boundaries_t(domain, &fluid->t)) {
      return 1;
    }
  }
  return 0;
}

