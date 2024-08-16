#include <string.h> // memset
#include "array.h"
#include "runge_kutta.h"

// coefficients of three-stage Runge-Kutta scheme | 8
static const double a0 = + 32. / 60., b0 =    0. / 60.;
static const double a1 = + 25. / 60., b1 = - 17. / 60.;
static const double a2 = + 45. / 60., b2 = - 25. / 60.;
const rkcoef_t rkcoefs[RKSTEPMAX] = {
  {.alpha = a0, .beta = b0, .gamma = a0 + b0},
  {.alpha = a1, .beta = b1, .gamma = a1 + b1},
  {.alpha = a2, .beta = b2, .gamma = a2 + b2},
};

int rkbuffers_init (
    const domain_t * const domain,
    const int nadds[NDIMS][2],
    rkbuffers_t * const buffers
) {
  if (0 != array.create(domain, nadds, sizeof(double), &buffers->alpha)) {
    return 1;
  }
  if (0 != array.create(domain, nadds, sizeof(double), &buffers->beta )) {
    return 1;
  }
  if (0 != array.create(domain, nadds, sizeof(double), &buffers->gamma)) {
    return 1;
  }
  return 0;
}

int rkbuffers_reset (
    rkbuffers_t * const buffers
) {
  array_t * const restrict alpha = &buffers->alpha;
  array_t * const restrict beta  = &buffers->beta ;
  array_t * const restrict gamma = &buffers->gamma;
  // stash previous RK source term,
  //   which is achieved by swapping
  //   the pointers to "data"
  double * const tmp = alpha->data;
  alpha->data = beta ->data;
  beta ->data = tmp;
  // zero-clear current RK source terms (exp/imp)
  // NOTE: when 0 == rkstep, rkcoef for "beta" is zero
  //   and thus zero-clearing"beta" buffers is not needed
  memset(alpha->data, 0, alpha->datasize);
  memset(gamma->data, 0, gamma->datasize);
  return 0;
}

