#include "runge_kutta.h"

const uint_fast8_t rk_a = 0;
const uint_fast8_t rk_b = 1;
const uint_fast8_t rk_g = 2;

// coefficients of three-stage Runge-Kutta scheme
static const double a0 = +32. / 60., b0 =   0. / 60.;
static const double a1 = +25. / 60., b1 = -17. / 60.;
static const double a2 = +45. / 60., b2 = -25. / 60.;
const rkcoef_t rkcoefs[RKSTEPMAX] = {
  {a0, b0, a0 + b0},
  {a1, b1, a1 + b1},
  {a2, b2, a2 + b2},
};

