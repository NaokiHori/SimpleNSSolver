#if !defined(FLUID_INTEGRATE_INTERNAL)
#define FLUID_INTEGRATE_INTERNAL

extern int compute_rhs_ux(
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_rhs_uy(
    const domain_t * domain,
    fluid_t * fluid
);

#if NDIMS == 3
extern int compute_rhs_uz(
    const domain_t * domain,
    fluid_t * fluid
);
#endif

extern int compute_rhs_t(
    const domain_t * domain,
    fluid_t * fluid
);

extern int buoyancy_ux(
    const domain_t * domain,
    fluid_t * fluid
);

extern int update_ux(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

extern int update_uy(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

#if NDIMS == 3
extern int update_uz(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);
#endif

extern int update_t(
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

#endif // FLUID_INTEGRATE_INTERNAL
