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
#include "array_macros/fluid/p.h"

static laplacians_t laplacians = {
  .is_initialised = false,
};

// map [1 : isize] in code to [0 : isize - 1] in memory
#define LAPX(I) lapx[(I) - 1]

static int init_laplacians(
    const domain_t * domain
) {
  // vector laplacian in x
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
  // vector laplacian in y
  {
    const double hy = domain->hy;
    const double l = 1. / hy / hy;
    const double u = 1. / hy / hy;
    const double c = - l - u;
    laplacians.lapy.l = l;
    laplacians.lapy.c = c;
    laplacians.lapy.u = u;
  }
  laplacians.is_initialised = true;
  return 0;
}

#define BEGIN \
  for (int cnt = 0, j = 1; j <= jsize; j++) { \
    for (int i = 1; i <= isize; i++, cnt++) {
#define END \
    } \
  }

// advected in x
static int advection_x (
    const domain_t * domain,
    const double * restrict uy,
    const double * restrict ux,
    double * restrict src
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict jdxc = domain->jdxc;
  BEGIN
    const double hx_xm = HXXF(i  );
    const double hx_xp = HXXF(i+1);
    const double jd_xm = JDXF(i  );
    const double jd_x0 = JDXC(i  );
    const double jd_xp = JDXF(i+1);
    const double ux_xm = + 0.5 * jd_xm / hx_xm * UX(i  , j-1)
                         + 0.5 * jd_xm / hx_xm * UX(i  , j  );
    const double ux_xp = + 0.5 * jd_xp / hx_xp * UX(i+1, j-1)
                         + 0.5 * jd_xp / hx_xp * UX(i+1, j  );
    const double l = - 0.5 * ux_xm;
    const double u = + 0.5 * ux_xp;
    const double c = - l - u;
    src[cnt] -= 1. / jd_x0 * (
        + l * UY(i-1, j  )
        + c * UY(i  , j  )
        + u * UY(i+1, j  )
    );
  END
  return 0;
}

// advected in y
static int advection_y (
    const domain_t * domain,
    const double * restrict uy,
    double * restrict src
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double hy = domain->hy;
  const double * restrict jdxc = domain->jdxc;
  BEGIN
    const double jd = JDXC(i  );
    const double uy_ym = + 0.5 * jd / hy * UY(i  , j-1)
                         + 0.5 * jd / hy * UY(i  , j  );
    const double uy_yp = + 0.5 * jd / hy * UY(i  , j  )
                         + 0.5 * jd / hy * UY(i  , j+1);
    const double l = - 0.5 * uy_ym;
    const double u = + 0.5 * uy_yp;
    const double c = - l - u;
    src[cnt] -= 1. / jd * (
        + l * UY(i  , j-1)
        + c * UY(i  , j  )
        + u * UY(i  , j+1)
    );
  END
  return 0;
}

// pressure gradient effect
static int pressure (
    const domain_t * domain,
    const double * restrict p,
    double * restrict src
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double hy = domain->hy;
  BEGIN
    src[cnt] -= 1. / hy * (
        - P(i  , j-1)
        + P(i  , j  )
    );
  END
  return 0;
}

// diffused in x
static int diffusion_x (
    const domain_t * domain,
    const double diffusivity,
    const double * restrict uy,
    double * restrict src
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const laplacian_t * lapx = laplacians.lapx;
  BEGIN
    src[cnt] += diffusivity * (
        + LAPX(i).l * UY(i-1, j  )
        + LAPX(i).c * UY(i  , j  )
        + LAPX(i).u * UY(i+1, j  )
    );
  END
  return 0;
}

// diffused in y
static int diffusion_y (
    const domain_t * domain,
    const double diffusivity,
    const double * restrict uy,
    double * restrict src
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const laplacian_t * lapy = &laplacians.lapy;
  BEGIN
    src[cnt] += diffusivity * (
        + (*lapy).l * UY(i  , j-1)
        + (*lapy).c * UY(i  , j  )
        + (*lapy).u * UY(i  , j+1)
    );
  END
  return 0;
}

// compute right-hand-side terms, which are added to buffers
int compute_rhs_uy(
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
  const double * restrict  p = fluid-> p.data;
  // buffer for explicit terms
  double * restrict srca = fluid->srcuy.alpha.data;
  // buffer for implicit terms
  double * restrict srcg = fluid->srcuy.gamma.data;
  const double diffusivity = fluid_compute_momentum_diffusivity(fluid);
  // advective contributions, always explicit
  advection_x(domain, uy, ux, srca);
  advection_y(domain, uy,     srca);
  // pressure-gradient contribution, always implicit
  pressure(domain, p, srcg);
  // diffusive contributions, can be explicit or implicit
  diffusion_x(domain, diffusivity, uy, param_m_implicit_x ? srcg : srca);
  diffusion_y(domain, diffusivity, uy, param_m_implicit_y ? srcg : srca);
  return 0;
}

// update y velocity field
int update_uy(
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
    };
    const size_t glsizes[NDIMS] = {
      domain->glsizes[0],
      domain->glsizes[1],
    };
    if (0 != linear_system_init(domain->info, implicit, glsizes, &linear_system)) {
      return 1;
    }
  }
  // compute increments
  {
    const double coef_a = rkcoefs[rkstep].alpha;
    const double coef_b = rkcoefs[rkstep].beta ;
    const double coef_g = rkcoefs[rkstep].gamma;
    const double * restrict srcuya = fluid->srcuy.alpha.data;
    const double * restrict srcuyb = fluid->srcuy.beta .data;
    const double * restrict srcuyg = fluid->srcuy.gamma.data;
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    double * restrict duy = linear_system.x1pncl;
    const size_t nitems = isize * jsize;
    for (size_t n = 0; n < nitems; n++) {
      duy[n] =
        + coef_a * dt * srcuya[n]
        + coef_b * dt * srcuyb[n]
        + coef_g * dt * srcuyg[n];
    }
  }
  // solve linear systems if necessary
  {
    const double prefactor =
      0.5 * rkcoefs[rkstep].gamma * dt * fluid_compute_momentum_diffusivity(fluid);
    // solve linear systems in x
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
    // solve linear systems in y
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
  }
  // update velocity field
  {
    const int isize = domain->mysizes[0];
    const int jsize = domain->mysizes[1];
    const double * restrict duy = linear_system.x1pncl;
    double * restrict uy = fluid->uy.data;
    BEGIN
      UY(i, j) += duy[cnt];
    END
    if (0 != fluid_update_boundaries_uy(domain, &fluid->uy)) {
      return 1;
    }
  }
  return 0;
}

