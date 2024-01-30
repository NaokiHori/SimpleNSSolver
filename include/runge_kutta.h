#if !defined(RUNGE_KUTTA_H)
#define RUNGE_KUTTA_H

#include <stdint.h>

// Runge-Kutta configurations
// indices
extern const uint_fast8_t rk_a; // 0
extern const uint_fast8_t rk_b; // 1
extern const uint_fast8_t rk_g; // 2
// NOTE: alpha, beta, gamma and thus three here
typedef double rkcoef_t[3];
// NOTE: only three-step Wray is allowed
#define RKSTEPMAX 3
extern const rkcoef_t rkcoefs[RKSTEPMAX];

#endif // RUNGE_KUTTA_H
