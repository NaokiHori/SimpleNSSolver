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
#include "internal.h"

// duxdx component
static int get_duxdx (
    const domain_t * domain,
    const fluid_t * fluid,
    double * quantity
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const double * restrict hxxc = domain->hxxc;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict ux = fluid->ux.data;
  const double diffusivity = fluid_compute_momentum_diffusivity(fluid);
  for (int j = 1; j <= jsize; j++) {
    for (int i = 1; i <= isize; i++) {
      const double hx = HXXC(i);
      const double jd = JDXC(i);
      const double dux =
        - UX(i  , j  )
        + UX(i+1, j  );
      *quantity += diffusivity * jd * pow(1. / hx * dux, 2.);
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
  const double hy = domain->hy;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict ux = fluid->ux.data;
  const double diffusivity = fluid_compute_momentum_diffusivity(fluid);
  for (int j = 1; j <= jsize; j++) {
    for (int i = 1; i <= isize + 1; i++) {
      const double jd = JDXF(i);
      const double dux =
        - UX(i  , j-1)
        + UX(i  , j  );
      *quantity += diffusivity * jd * pow(1. / hy * dux, 2.);
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
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  const double * restrict uy = fluid->uy.data;
  const double diffusivity = fluid_compute_momentum_diffusivity(fluid);
  for (int j = 1; j <= jsize; j++) {
    for (int i = 1; i <= isize + 1; i++) {
      const double hx = HXXF(i);
      const double jd = JDXF(i);
      const double duy =
        - UY(i-1, j  )
        + UY(i  , j  );
      *quantity += diffusivity * jd * pow(1. / hx * duy, 2.);
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
  const double hy = domain->hy;
  const double * restrict jdxc = domain->jdxc;
  const double * restrict uy = fluid->uy.data;
  const double diffusivity = fluid_compute_momentum_diffusivity(fluid);
  for (int j = 1; j <= jsize; j++) {
    for (int i = 1; i <= isize; i++) {
      const double jd = JDXC(i);
      const double duy =
        - UY(i  , j  )
        + UY(i  , j+1);
      *quantity += diffusivity * jd * pow(1. / hy * duy, 2.);
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
  get_duydx(domain, fluid, &quantity);
  get_duydy(domain, fluid, &quantity);
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

