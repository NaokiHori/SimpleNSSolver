#if !defined(PARAM_H)
#define PARAM_H

// fixed parameters, which are usually fixed
//   but still user can easily control, are declared
// they are defined under src/param/xxx.c

#include <stdbool.h>

/* buoyancy.c */
// flag to specify whether the buoyancy force is added
//   to the RHS of the momentum equation or not
// if not, the temperature behaves as a passive scalar
extern const bool param_add_buoyancy;

/* implicit.c */
// flags to specify the diffusive treatment of the momentum equations
extern const bool param_m_implicit_x;
extern const bool param_m_implicit_y;
// flags to specify the diffusive treatment of the temperature equation
extern const bool param_t_implicit_x;
extern const bool param_t_implicit_y;

/* boundary-condition.c */
// NOTE: changing values may break the Nusselt balance
// NOTE: impermeable walls and Neumann BC for the pressure are unchangeable
// negative-x-wall velocity in y direction
extern const double param_uy_xm;
// positive-x-wall velocity in y direction
extern const double param_uy_xp;
// negative-x-wall temperature
extern const double param_t_xm;
// positive-x-wall temperature
extern const double param_t_xp;

#endif // PARAM_H
