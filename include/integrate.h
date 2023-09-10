#if !defined(INTEGRATE_H)
#define INTEGRATE_H

// main integrator
extern int integrate(
    const domain_t * domain,
    fluid_t * fluid,
    double * dt
);

#endif // INTEGRATE_H
