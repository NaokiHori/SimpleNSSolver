#if !defined(LOGGING_NUSSELT_INTERNAL_H)
#define LOGGING_NUSSELT_INTERNAL_H

extern double logging_internal_compute_reference_heat_flux(
    const domain_t * domain,
    const fluid_t * fluid
);

extern double logging_internal_compute_nu_heat_flux(
    const domain_t * domain,
    const fluid_t * fluid
);

extern double logging_internal_compute_nu_kinetic_energy_injection(
    const domain_t * domain,
    const fluid_t * fluid
);

extern double logging_internal_compute_nu_kinetic_energy_dissipation(
    const domain_t * domain,
    const fluid_t * fluid
);

extern double logging_internal_compute_nu_thermal_energy_dissipation(
    const domain_t * domain,
    const fluid_t * fluid
);

#endif // LOGGING_NUSSELT_INTERNAL_H
