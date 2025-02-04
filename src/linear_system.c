#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "memory.h"
#include "sdecomp.h"
#include "domain.h"
#include "linear_system.h"
#include "tdm.h"

static size_t get_nitems(
    const size_t sizes[NDIMS]
){
  // helper function just to multiply sizes and return
  size_t nitems = 1;
  for(int dim = 0; dim < NDIMS; dim++){
    nitems *= sizes[dim];
  }
  return nitems;
}

/**
 * @brief initialise linear solver to update field implicitly
 * @param[in] info     : information about domain decomposition
 * @param[in] implicit : treatment of the diffusive terms in each direction
 * @param[in] glsizes  : GLOBAL size of array
 * @return             : structure storing buffers and plans to solve linear systems in each direction
 */
int linear_system_init(
    const sdecomp_info_t * info,
    const bool implicit[NDIMS],
    const size_t glsizes[NDIMS],
    linear_system_t * linear_system
){
  if(linear_system->is_initialised){
    printf("this linear_system object is already initialised\n");
    return 1;
  }
  memcpy(linear_system->implicit, implicit, sizeof(bool) * NDIMS);
  // pencils (and their sizes) to store input and output of linear systems
  double * restrict * x1pncl = &linear_system->x1pncl;
  double * restrict * y1pncl = &linear_system->y1pncl;
  double * restrict * z2pncl = &linear_system->z2pncl;
  // check size first
  size_t * x1pncl_mysizes = linear_system->x1pncl_mysizes;
  size_t * y1pncl_mysizes = linear_system->y1pncl_mysizes;
  size_t * z2pncl_mysizes = linear_system->z2pncl_mysizes;
  for(int dim = 0; dim < NDIMS; dim++){
    if(0 != sdecomp.get_pencil_mysize(info, SDECOMP_X1PENCIL, dim, glsizes[dim], x1pncl_mysizes + dim)) return 1;
    if(0 != sdecomp.get_pencil_mysize(info, SDECOMP_Y1PENCIL, dim, glsizes[dim], y1pncl_mysizes + dim)) return 1;
    if(0 != sdecomp.get_pencil_mysize(info, SDECOMP_Z2PENCIL, dim, glsizes[dim], z2pncl_mysizes + dim)) return 1;
  }
  // allocate pencils if needed
  // NOTE: x1pncl is not needed for fully-explicit case,
  //   but I always allocate it here for simplicity (to store delta values)
  *x1pncl =               memory_calloc(get_nitems(x1pncl_mysizes), sizeof(double));
  *y1pncl = implicit[1] ? memory_calloc(get_nitems(y1pncl_mysizes), sizeof(double)) : NULL;
  *z2pncl = implicit[2] ? memory_calloc(get_nitems(z2pncl_mysizes), sizeof(double)) : NULL;
  // initialise parallel matrix transpose if needed
  // between x1 and y1
  if(implicit[1]){
    sdecomp_transpose_plan_t ** plan_f = &linear_system->transposer_x1_to_y1;
    sdecomp_transpose_plan_t ** plan_b = &linear_system->transposer_y1_to_x1;
    if(0 != sdecomp.transpose.construct(info, SDECOMP_X1PENCIL, SDECOMP_Y1PENCIL, glsizes, sizeof(double), plan_f)) return 1;
    if(0 != sdecomp.transpose.construct(info, SDECOMP_Y1PENCIL, SDECOMP_X1PENCIL, glsizes, sizeof(double), plan_b)) return 1;
  }
  // between x1 and z2
  if(implicit[2]){
    sdecomp_transpose_plan_t ** plan_f = &linear_system->transposer_x1_to_z2;
    sdecomp_transpose_plan_t ** plan_b = &linear_system->transposer_z2_to_x1;
    if(0 != sdecomp.transpose.construct(info, SDECOMP_X1PENCIL, SDECOMP_Z2PENCIL, glsizes, sizeof(double), plan_f)) return 1;
    if(0 != sdecomp.transpose.construct(info, SDECOMP_Z2PENCIL, SDECOMP_X1PENCIL, glsizes, sizeof(double), plan_b)) return 1;
  }
  // initialise tri-diagonal matrix solvers
  // Thomas algorithm in x direction
  if(implicit[0]){
    if(0 != tdm.construct(
        /* size of system */ (int)(x1pncl_mysizes[0]),
        /* number of rhs  */ (int)(x1pncl_mysizes[1] * x1pncl_mysizes[2]),
        /* is periodic    */ false,
        /* is complex     */ false,
        /* output         */ &linear_system->tdm_x
    )) return 1;
  }
  // Thomas algorithm in y direction
  if(implicit[1]){
    if(0 != tdm.construct(
        /* size of system */ (int)(y1pncl_mysizes[1]),
        /* number of rhs  */ (int)(y1pncl_mysizes[2] * y1pncl_mysizes[0]),
        /* is periodic    */ true,
        /* is complex     */ false,
        /* output         */ &linear_system->tdm_y
    )) return 1;
  }
  // Thomas algorithm in z direction
  if(implicit[2]){
    if(0 != tdm.construct(
        /* size of system */ (int)(z2pncl_mysizes[2]),
        /* number of rhs  */ (int)(z2pncl_mysizes[0] * z2pncl_mysizes[1]),
        /* is periodic    */ true,
        /* is complex     */ false,
        /* output         */ &linear_system->tdm_z
    )) return 1;
  }
  linear_system->is_initialised = true;
  return 0;
}

int linear_system_finalise(
    linear_system_t * linear_system
){
  if(NULL == linear_system || !linear_system->is_initialised){
    return 1;
  }
  const bool * implicit = linear_system->implicit;
  memory_free(linear_system->x1pncl);
  if(implicit[0]){
    tdm.destruct(linear_system->tdm_x);
  }
  if(implicit[1]){
    memory_free(linear_system->y1pncl);
    sdecomp.transpose.destruct(linear_system->transposer_x1_to_y1);
    sdecomp.transpose.destruct(linear_system->transposer_y1_to_x1);
    tdm.destruct(linear_system->tdm_y);
  }
  if(implicit[2]){
    memory_free(linear_system->z2pncl);
    sdecomp.transpose.destruct(linear_system->transposer_x1_to_z2);
    sdecomp.transpose.destruct(linear_system->transposer_z2_to_x1);
    tdm.destruct(linear_system->tdm_z);
  }
  return 0;
}

