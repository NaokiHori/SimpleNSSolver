#if !defined(DECIDE_DT_H)
#define DECIDE_DT_H

// decide next time step size
extern int decide_dt(
    const domain_t * domain,
    const fluid_t * fluid,
    double * dt
);

#endif // DECIDE_DT_H
