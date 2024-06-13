#if !defined(LOGGING_INTERNAL_H)
#define LOGGING_INTERNAL_H

extern int logging_internal_output (
    const char fname[],
    const domain_t * domain,
    const double time,
    const size_t nitems,
    const double * items
);

extern int logging_check_max_divergence (
    const char fname[],
    const domain_t * domain,
    const double time,
    const fluid_t * fluid
);

extern int logging_check_total_momentum (
    const char fname[],
    const domain_t * domain,
    const double time,
    const fluid_t * fluid
);

extern int logging_check_total_energy (
    const char fname[],
    const domain_t * domain,
    const double time,
    const fluid_t * fluid
);

extern int logging_check_heat_transfer (
    const char fname[],
    const domain_t * domain,
    const double time,
    const fluid_t * fluid
);

extern int logging_check_injected_squared_velocity (
    const char fname[],
    const domain_t * domain,
    const double time,
    const fluid_t * fluid
);

int logging_check_dissipated_squared_velocity (
    const char fname[],
    const domain_t * domain,
    const double time,
    const fluid_t * fluid
);

int logging_check_dissipated_squared_temperature (
    const char fname[],
    const domain_t * domain,
    const double time,
    const fluid_t * fluid
);

#endif // LOGGING_INTERNAL_H
