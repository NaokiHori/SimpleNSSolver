#if !defined(FLUID_PREDICT_INTERNAL)
#define FLUID_PREDICT_INTERNAL

#include <stdbool.h>

// store approximation of laplacian
typedef struct {
  double l;
  double c;
  double u;
} laplacian_t;

// store Laplacian for each directoin
typedef struct {
  bool is_initialised;
  laplacian_t * lapx;
  laplacian_t lapy;
#if NDIMS == 3
  laplacian_t lapz;
#endif
} laplacians_t;

extern int compute_rhs_ux (
    const domain_t * domain,
    fluid_t * fluid
);

extern int compute_rhs_uy (
    const domain_t * domain,
    fluid_t * fluid
);

#if NDIMS == 3
extern int compute_rhs_uz (
    const domain_t * domain,
    fluid_t * fluid
);
#endif

extern int compute_rhs_t (
    const domain_t * domain,
    fluid_t * fluid
);

extern int update_ux (
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

extern int update_uy (
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

#if NDIMS == 3
extern int update_uz (
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);
#endif

extern int update_t (
    const domain_t * domain,
    const size_t rkstep,
    const double dt,
    fluid_t * fluid
);

#endif // FLUID_PREDICT_INTERNAL
