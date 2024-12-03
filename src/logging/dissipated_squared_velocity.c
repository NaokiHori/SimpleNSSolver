#include <math.h>
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/hxxc.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/domain/jdxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/uz.h"
#include "internal.h"

// duxdx component
static int get_duxdx (
    const domain_t * domain,
    const fluid_t * fluid,
    double * quantity
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hxxc = domain->hxxc;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict ux = fluid->ux.data;
  const double diffusivity = fluid_compute_momentum_diffusivity(fluid);
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize; i++) {
        const double hx = HXXC(i);
        const double jd = JDXC(i);
        const double dux =
          - UX(i  , j  , k  )
          + UX(i+1, j  , k  );
        *quantity += diffusivity * jd * pow(1. / hx * dux, 2.);
      }
    }
  }
  return 0;
}

// duxdy component
static int get_duxdy (
    const domain_t * domain,
    const fluid_t * fluid,
    double * quantity
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hy = domain->hy;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict ux = fluid->ux.data;
  const double diffusivity = fluid_compute_momentum_diffusivity(fluid);
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize + 1; i++) {
        const double jd = JDXF(i);
        const double dux =
          - UX(i  , j-1, k  )
          + UX(i  , j  , k  );
        *quantity += diffusivity * jd * pow(1. / hy * dux, 2.);
      }
    }
  }
  return 0;
}

// duxdz component
static int get_duxdz (
    const domain_t * domain,
    const fluid_t * fluid,
    double * quantity
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict ux = fluid->ux.data;
  const double diffusivity = fluid_compute_momentum_diffusivity(fluid);
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize + 1; i++) {
        const double jd = JDXF(i);
        const double dux =
          - UX(i  , j  , k-1)
          + UX(i  , j  , k  );
        *quantity += diffusivity * jd * pow(1. / hz * dux, 2.);
      }
    }
  }
  return 0;
}

// duydx component
static int get_duydx (
    const domain_t * domain,
    const fluid_t * fluid,
    double * quantity
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict uy = fluid->uy.data;
  const double diffusivity = fluid_compute_momentum_diffusivity(fluid);
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize + 1; i++) {
        const double hx = HXXF(i);
        const double jd = JDXF(i);
        const double duy =
          - UY(i-1, j  , k  )
          + UY(i  , j  , k  );
        *quantity += diffusivity * jd * pow(1. / hx * duy, 2.);
      }
    }
  }
  return 0;
}

// duydy component
static int get_duydy (
    const domain_t * domain,
    const fluid_t * fluid,
    double * quantity
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hy = domain->hy;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict uy = fluid->uy.data;
  const double diffusivity = fluid_compute_momentum_diffusivity(fluid);
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize; i++) {
        const double jd = JDXC(i);
        const double duy =
          - UY(i  , j  , k  )
          + UY(i  , j+1, k  );
        *quantity += diffusivity * jd * pow(1. / hy * duy, 2.);
      }
    }
  }
  return 0;
}

// duydz component
static int get_duydz (
    const domain_t * domain,
    const fluid_t * fluid,
    double * quantity
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict uy = fluid->uy.data;
  const double diffusivity = fluid_compute_momentum_diffusivity(fluid);
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize; i++) {
        const double jd = JDXC(i);
        const double duy =
          - UY(i  , j  , k-1)
          + UY(i  , j  , k  );
        *quantity += diffusivity * jd * pow(1. / hz * duy, 2.);
      }
    }
  }
  return 0;
}

// duzdx component
static int get_duzdx (
    const domain_t * domain,
    const fluid_t * fluid,
    double * quantity
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict uz = fluid->uz.data;
  const double diffusivity = fluid_compute_momentum_diffusivity(fluid);
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize + 1; i++) {
        const double hx = HXXF(i);
        const double jd = JDXF(i);
        const double duz =
          - UZ(i-1, j  , k  )
          + UZ(i  , j  , k  );
        *quantity += diffusivity * jd * pow(1. / hx * duz, 2.);
      }
    }
  }
  return 0;
}

// duzdy component
static int get_duzdy (
    const domain_t * domain,
    const fluid_t * fluid,
    double * quantity
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hy = domain->hy;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict uz = fluid->uz.data;
  const double diffusivity = fluid_compute_momentum_diffusivity(fluid);
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize; i++) {
        const double jd = JDXC(i);
        const double duz =
          - UZ(i  , j-1, k  )
          + UZ(i  , j  , k  );
        *quantity += diffusivity * jd * pow(1. / hy * duz, 2.);
      }
    }
  }
  return 0;
}

// duzdz component
static int get_duzdz (
    const domain_t * domain,
    const fluid_t * fluid,
    double * quantity
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double hz = domain->hz;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict uz = fluid->uz.data;
  const double diffusivity = fluid_compute_momentum_diffusivity(fluid);
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize; i++) {
        const double jd = JDXC(i);
        const double duz =
          - UZ(i  , j  , k  )
          + UZ(i  , j  , k+1);
        *quantity += diffusivity * jd * pow(1. / hz * duz, 2.);
      }
    }
  }
  return 0;
}

/**
 * @brief compute dissipation of squared velocity
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] time   : current simulation time
 * @param[in] fluid  : diffusivity and velocity field
 * @return           : error code
 */
int logging_check_dissipated_squared_velocity (
    const char fname[],
    const domain_t * domain,
    const double time,
    const fluid_t * fluid
) {
  const int root = 0;
  int myrank = root;
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_rank(domain->info, &myrank);
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  double quantity = 0.;
  get_duxdx(domain, fluid, &quantity);
  get_duxdy(domain, fluid, &quantity);
  get_duxdz(domain, fluid, &quantity);
  get_duydx(domain, fluid, &quantity);
  get_duydy(domain, fluid, &quantity);
  get_duydz(domain, fluid, &quantity);
  get_duzdx(domain, fluid, &quantity);
  get_duzdy(domain, fluid, &quantity);
  get_duzdz(domain, fluid, &quantity);
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : &quantity;
  void * recvbuf = &quantity;
  MPI_Reduce(sendbuf, recvbuf, 1, MPI_DOUBLE, MPI_SUM, root, comm_cart);
  logging_internal_output(
      fname,
      domain,
      time,
      1,
      &quantity
  );
  return 0;
}

