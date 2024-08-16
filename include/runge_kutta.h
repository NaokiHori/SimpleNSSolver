#if !defined(RUNGE_KUTTA_H)
#define RUNGE_KUTTA_H

#include "domain.h"
#include "array.h"

// Runge-Kutta configurations
// NOTE: only three-step Wray is allowed
#define RKSTEPMAX 3

// coefficients in front of dt
typedef struct {
  double alpha;
  double beta;
  double gamma;
} rkcoef_t;
extern const rkcoef_t rkcoefs[RKSTEPMAX];

// buffers to store
//   explicit terms (current: alpha, previous: beta)
//   and
//   implicit terms (gamma)
typedef struct {
  array_t alpha;
  array_t beta;
  array_t gamma;
} rkbuffers_t;

extern int rkbuffers_init (
    const domain_t * const domain,
    const int nadds[NDIMS][2],
    rkbuffers_t * const buffers
);

extern int rkbuffers_reset (
    rkbuffers_t * const buffers
);

#endif // RUNGE_KUTTA_H
