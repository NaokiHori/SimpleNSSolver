#if !defined(DOMAIN_H)
#define DOMAIN_H

#include "sdecomp.h"

// definition of a structure domain_t | 25
/**
 * @struct domain_t
 * @brief struct storing parameters relevant to spatial domain
 * @var info       : MPI domain decomposition
 * @var glsizes    : global     number of grid points in each direction
 * @var mysizes    : local (my) number of grid points in each direction
 * @var offsets    : offsets to my starting index in each direction
 * @var lengths    : domain size in each direction
 * @var xf, xc     : cell-face and cell-center locations in x direction
 * @var hxxf, hxxc : wall-normal scale factors at faces and centers
 * @var hy, hz     : scale factors in the homogeneous directions
 */
typedef struct {
  sdecomp_info_t * info;
  size_t glsizes[NDIMS];
  size_t mysizes[NDIMS];
  size_t offsets[NDIMS];
  double lengths[NDIMS];
  double * restrict xf, * restrict xc;
  double * restrict hxxf, * restrict hxxc;
  double hy;
#if NDIMS == 3
  double hz;
#endif
  double * restrict jdxf, * restrict jdxc;
} domain_t;

// constructor
extern int domain_init(
    const char dirname_ic[],
    domain_t * domain
);

// save members which are necessary to restart
extern int domain_save(
    const char dirname[],
    const domain_t * domain
);

#endif // DOMAIN_H
