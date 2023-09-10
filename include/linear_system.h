#if !defined(LINEAR_SYSTEM_H)
#define LINEAR_SYSTEM_H

#include <stdbool.h>
#include "sdecomp.h"
#include "domain.h"
#include "tdm.h"

/**
 * @struct linear_system_t
 * @brief structure storing buffers and plans to solve tri-diagonal linear systems in each dimension A x = b
 * @var is_initialised      : flag to check the variable is initialised
 * @var implicit            : flags whether directions are treated implicitly,
 *                              namely linear systems are to be solved
 * @var x1pncl              : buffers to store x1-pencil
 * @var y1pncl              : buffers to store y1-pencil
 * @var z2pncl              : buffers to store z2-pencil
 * @var x1pncl_mysizes      : size of (local) x1pencil
 * @var y1pncl_mysizes      : size of (local) y1pencil
 * @var z2pncl_mysizes      : size of (local) z2pencil
 * @var tdm_[x-z]           : thomas algorithm solvers in all directions
 * @var transposer_xx_to_xx : plans to transpose between two pencils
 */
typedef struct {
  bool is_initialised;
  bool implicit[NDIMS];
  double * restrict x1pncl;
  double * restrict y1pncl;
#if NDIMS == 3
  double * restrict z2pncl;
#endif
  size_t x1pncl_mysizes[NDIMS];
  size_t y1pncl_mysizes[NDIMS];
#if NDIMS == 3
  size_t z2pncl_mysizes[NDIMS];
#endif
  tdm_info_t * tdm_x;
  tdm_info_t * tdm_y;
#if NDIMS == 3
  tdm_info_t * tdm_z;
#endif
  sdecomp_transpose_plan_t * transposer_x1_to_y1;
  sdecomp_transpose_plan_t * transposer_y1_to_x1;
#if NDIMS == 3
  sdecomp_transpose_plan_t * transposer_x1_to_z2;
  sdecomp_transpose_plan_t * transposer_z2_to_x1;
#endif
} linear_system_t;

extern int linear_system_init(
    const sdecomp_info_t * info,
    const bool implicit[NDIMS],
    const size_t glsizes[NDIMS],
    linear_system_t * linear_system
);

extern int linear_system_finalise(
    linear_system_t * linear_system
);

#endif // LINEAR_SYSTEM_H
