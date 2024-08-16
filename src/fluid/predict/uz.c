#if NDIMS == 3
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
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/uz.h"
#include "array_macros/fluid/p.h"

static laplacians_t laplacians = {
  .is_initialised = false,
};

// map [1 : isize] in code to [0 : isize - 1] in memory
#define LAPX(I) lapx[(I) - 1]

static int init_laplacians (
    const domain_t * domain
) {
  // vector laplacian in x | 15
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
  // vector laplacian in y | 9
  {
    const double hy = domain->hy;
    const double l = 1. / hy / hy;
    const double u = 1. / hy / hy;
    const double c = - l - u;
    laplacians.lapy.l = l;
    laplacians.lapy.c = c;
    laplacians.lapy.u = u;
  }
  // vector laplacian in z | 9
  {
    const double hz = domain->hz;
    const double l = 1. / hz / hz;
    const double u = 1. / hz / hz;
    const double c = - l - u;
    laplacians.lapz.l = l;
    laplacians.lapz.c = c;
    laplacians.lapz.u = u;
  }
  laplacians.is_initialised = true;
  return 0;
}

#define BEGIN \
  for (int cnt = 0, k = 1; k <= ksize; k++) { \
    for (int j = 1; j <= jsize; j++) { \
      for (int i = 1; i <= isize; i++, cnt++) {
#define END \
      } \
    } \
  }

// advected in x | 33
static int advection_x (
    const domain_t * domain,
    const double * restrict uz,
    const double * restrict ux,
    double * restrict src
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  BEGIN
    const double hx_xm = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i  );
    const double jd_x0 = JDXC(i  );
    const double jd_xp = JDXF(i+1);
    const double ux_xm = + 0.5 * jd_xm / hx_xm * UX(i  , j  , k-1)
                         + 0.5 * jd_xm / hx_xm * UX(i  , j  , k  );
    const double ux_xp = + 0.5 * jd_xp / hx_xp * UX(i+1, j  , k-1)
                         + 0.5 * jd_xp / hx_xp * UX(i+1, j  , k  );
    const double l = - 0.5 * ux_xm;
    const double u = + 0.5 * ux_xp;
    const double c = - l - u;
    src[cnt] -= 1. / jd_x0 * (
        + l * UZ(i-1, j  , k  )
        + c * UZ(i  , j  , k  )
        + u * UZ(i+1, j  , k  )
    );
  END
  return 0;
}

// advected in y | 28
static int advection_y (
    const domain_t * domain,
    const double * restrict uz,
    const double * restrict uy,
    double * restrict src
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hy = domain->hy;
  const double * restrict jdxc = domain->jdxc;
  BEGIN
    const double jd = JDXC(i  );
    const double uy_ym = + 0.5 * jd / hy * UY(i  , j  , k-1)
                         + 0.5 * jd / hy * UY(i  , j  , k  );
    const double uy_yp = + 0.5 * jd / hy * UY(i  , j+1, k-1)
                         + 0.5 * jd / hy * UY(i  , j+1, k  );
    const double l = - 0.5 * uy_ym;
    const double u = + 0.5 * uy_yp;
    const double c = - l - u;
    src[cnt] -= 1. / jd * (
        + l * UZ(i  , j-1, k  )
        + c * UZ(i  , j  , k  )
        + u * UZ(i  , j+1, k  )
    );
  END
  return 0;
}

// advected in z | 27
static int advection_z (
    const domain_t * domain,
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
    const double uz_zm = + 0.5 * jd / hz * UZ(i  , j  , k-1)
                         + 0.5 * jd / hz * UZ(i  , j  , k  );
    const double uz_zp = + 0.5 * jd / hz * UZ(i  , j  , k  )
                         + 0.5 * jd / hz * UZ(i  , j  , k+1);
    const double l = - 0.5 * uz_zm;
    const double u = + 0.5 * uz_zp;
    const double c = - l - u;
    src[cnt] -= 1. / jd * (
        + l * UZ(i  , j  , k-1)
        + c * UZ(i  , j  , k  )
        + u * UZ(i  , j  , k+1)
    );
  END
  return 0;
}

// pressure gradient effect | 17
static int pressure (
    const domain_t * domain,
    const double * restrict p,
    double * restrict src
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  BEGIN
    src[cnt] -= 1. / hz * (
        - P(i  , j  , k-1)
        + P(i  , j  , k  )
    );
  END
  return 0;
}

// diffused in x | 19
static int diffusion_x (
    const domain_t * domain,
    const double diffusivity,
    const double * restrict uz,
    double * restrict src
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t * lapx = laplacians.lapx;
  BEGIN
    src[cnt] += diffusivity * (
        + LAPX(i).l * UZ(i-1, j  , k  )
        + LAPX(i).c * UZ(i  , j  , k  )
        + LAPX(i).u * UZ(i+1, j  , k  )
    );
  END
  return 0;
}

// diffused in y | 19
static int diffusion_y (
    const domain_t * domain,
    const double diffusivity,
    const double * restrict uz,
    double * restrict src
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t * lapy = &laplacians.lapy;
  BEGIN
    src[cnt] += diffusivity * (
        + (*lapy).l * UZ(i  , j-1, k  )
        + (*lapy).c * UZ(i  , j  , k  )
        + (*lapy).u * UZ(i  , j+1, k  )
    );
  END
  return 0;
}

// diffused in z | 19
static int diffusion_z (
    const domain_t * domain,
    const double diffusivity,
    const double * restrict uz,
    double * restrict src
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const laplacian_t * lapz = &laplacians.lapz;
  BEGIN
    src[cnt] += diffusivity * (
        + (*lapz).l * UZ(i  , j  , k-1)
        + (*lapz).c * UZ(i  , j  , k  )
        + (*lapz).u * UZ(i  , j  , k+1)
    );
  END
  return 0;
}

// compute right-hand-side terms, which are added to buffers | 30
int compute_rhs_uz (
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
  const double * restrict uz = fluid->uz.data;
  const double * restrict  p = fluid-> p.data;
  // buffer for explicit terms
  double * restrict srca = fluid->srcuz.alpha.data;
  // buffer for implicit terms
  double * restrict srcg = fluid->srcuz.gamma.data;
  const double diffusivity = fluid_compute_momentum_diffusivity(fluid);
  // advective contributions, always explicit
  advection_x(domain, uz, ux, srca);
  advection_y(domain, uz, uy, srca);
  advection_z(domain, uz,     srca);
  // pressure-gradient contribution, always implicit
  pressure(domain, p, srcg);
  // diffusive contributions, can be explicit or implicit
  diffusion_x(domain, diffusivity, uz, param_m_implicit_x ? srcg : srca);
  diffusion_y(domain, diffusivity, uz, param_m_implicit_y ? srcg : srca);
  diffusion_z(domain, diffusivity, uz, param_m_implicit_z ? srcg : srca);
  return 0;
}

// update z velocity field
int update_uz (
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
    //   for implicit diffusive term treatment
    const bool implicit[NDIMS] = {
      param_m_implicit_x,
      param_m_implicit_y,
      param_m_implicit_z,
    };
    const size_t glsizes[NDIMS] = {
      domain->glsizes[0],
      domain->glsizes[1],
      domain->glsizes[2],
    };
    if (0 != linear_system_init(domain->info, implicit, glsizes, &linear_system)) {
      return 1;
    }
  }
  // compute increments | 19
  {
    const double coef_a = rkcoefs[rkstep].alpha;
    const double coef_b = rkcoefs[rkstep].beta ;
    const double coef_g = rkcoefs[rkstep].gamma;
    const double * restrict srcuza = fluid->srcuz.alpha.data;
    const double * restrict srcuzb = fluid->srcuz.beta .data;
    const double * restrict srcuzg = fluid->srcuz.gamma.data;
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const int ksize = domain->mysizes[2];
    double * restrict duz = linear_system.x1pncl;
    const size_t nitems = isize * jsize * ksize;
    for (size_t n = 0; n < nitems; n++) {
      duz[n] =
        + coef_a * dt * srcuza[n]
        + coef_b * dt * srcuzb[n]
        + coef_g * dt * srcuzg[n];
    }
  }
  // solve linear systems if necessary
  {
    const double prefactor =
      0.5 * rkcoefs[rkstep].gamma * dt * fluid_compute_momentum_diffusivity(fluid);
    // solve linear systems in x | 18
    if (param_m_implicit_x) {
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
      tdm.solve(tdm_info, linear_system.x1pncl);
    }
    // solve linear systems in y | 28
    if (param_m_implicit_y) {
      sdecomp.transpose.execute(
          linear_system.transposer_x1_to_y1,
          linear_system.x1pncl,
          linear_system.y1pncl
      );
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
      tdm.solve(tdm_info, linear_system.y1pncl);
      sdecomp.transpose.execute(
          linear_system.transposer_y1_to_x1,
          linear_system.y1pncl,
          linear_system.x1pncl
      );
    }
    // solve linear systems in z | 28
    if (param_m_implicit_z) {
      sdecomp.transpose.execute(
          linear_system.transposer_x1_to_z2,
          linear_system.x1pncl,
          linear_system.z2pncl
      );
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
      tdm.solve(tdm_info, linear_system.z2pncl);
      sdecomp.transpose.execute(
          linear_system.transposer_z2_to_x1,
          linear_system.z2pncl,
          linear_system.x1pncl
      );
    }
  }
  // update velocity field | 13
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const int ksize = domain->mysizes[2];
    const double * restrict duz = linear_system.x1pncl;
    double * restrict uz = fluid->uz.data;
    BEGIN
      UZ(i, j, k) += duz[cnt];
    END
    if (0 != fluid_update_boundaries_uz(domain, &fluid->uz)) {
      return 1;
    }
  }
  return 0;
}
#endif
