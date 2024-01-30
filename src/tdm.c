#include <math.h>
#include <float.h>
#include <complex.h>
#include <fftw3.h>
#include "memory.h"
#include "tdm.h"

/**
 * @brief kernel function to solve a linear system
 * @param[in]    n : matrix size
 * @param[in]    l : lower  diagonal part
 * @param[in]    c : center diagonal part
 * @param[in]    u : upper  diagonal part
 * @param[inout] q : right-hand-side & answers
 * @return         : error code
 */
#define GTSV(type) \
  static int gtsv_##type( \
      const int n, \
      const double * restrict l, \
      const double * restrict c, \
      const double * restrict u, \
      double * restrict v, \
      type * restrict q \
){ \
    /* divide the first row by center-diagonal term | 2 */ \
    v[0] = u[0] / c[0]; \
    q[0] = q[0] / c[0]; \
    /* forward substitution | 7 */ \
    for(int i = 1; i < n - 1; i++){ \
      /* assume positive-definite system */ \
      /*   to skip zero-division checks */ \
      double val = 1. / (c[i] - l[i] * v[i-1]); \
      v[i] = val *      (u[i]                ); \
      q[i] = val *      (q[i] - l[i] * q[i-1]); \
    } \
    /* last row, do the same thing but consider singularity | 7 */ \
    double val = c[n-1] - l[n-1] * v[n-2]; \
    if(fabs(val) > DBL_EPSILON){ \
      q[n-1] = 1. / val * (q[n-1] - l[n-1] * q[n-2]); \
    }else{ \
      /* singular, zero mean */ \
      q[n-1] = 0.; \
    } \
    /* backward substitution | 3 */ \
    for(int i = n - 2; i >= 0; i--){ \
      q[i] -= v[i] * q[i+1]; \
    } \
    return 0; \
  }

/**
 * @brief solve linear system
 * @param[in]    n           : size of tri-diagonal matrix
 * @param[in]    m           : how many right-hand-sides do you want to solve?
 * @param[in]    is_periodic : periodic boundary condition is imposed
 *                               (Sherman-Morrison formula is used)
 *                               or not (normal Thomas algorithm is used)
 * @param[in]    l           : pointer to lower- diagonal components
 * @param[in]    c           : pointer to center-diagonal components
 * @param[in]    u           : pointer to upper- diagonal components
 * @param[inout] q           : right-hand-sides (size: "n", repeat for "m" times) & answers
 *                               N.B. memory is contiguous in "n" direction, sparse in "m" direction
 * @return                   : error code
 */
#define TDM_SOLVE(type) \
  static int tdm_solve_##type( \
      const int n, \
      const int nrhs, \
      const bool is_periodic, \
      const double * restrict l, \
      const double * restrict c, \
      const double * restrict u, \
      double * restrict v, \
      type * restrict q, \
      double * restrict q1 \
){ \
    if(is_periodic){ \
      /* solve additional system coming from periodicity | 7 */ \
      for(int i = 0; i < n-1; i++){ \
        q1[i] \
          = i ==   0 ? -l[i] \
          : i == n-2 ? -u[i] \
          : 0.; \
      } \
      gtsv_double(n-1, l, c, u, v, q1); \
      for(int j = 0; j < nrhs; j++){ \
        /* solve normal system | 2 */ \
        type *q0 = q + j * n; \
        gtsv_##type(n-1, l, c, u, v, q0); \
        /* find x_{n-1} | 3 */ \
        type   num = q0[n-1] - u[n-1] * q0[0] - l[n-1] * q0[n-2]; \
        double den = c [n-1] + u[n-1] * q1[0] + l[n-1] * q1[n-2]; \
        q0[n-1] = fabs(den) < DBL_EPSILON ? 0. : num / den; \
        /* solve original system | 3 */ \
        for(int i = 0; i < n-1; i++){ \
          q0[i] = q0[i] + q0[n-1] * q1[i]; \
        } \
      } \
    }else{ \
      for(int j = 0; j < nrhs; j++){ \
        gtsv_##type(n, l, c, u, v, q + j * n); \
      } \
    } \
    return 0; \
  }

// expand macros to define solvers
GTSV(double)
GTSV(fftw_complex)
TDM_SOLVE(double)
TDM_SOLVE(fftw_complex)

// definition of tdm_info_t_
/**
 * @struct tdm_info_t_
 * @brief struct to keep information about tri-diagonal system and internal buffers
 * @var size         : size of the system
 * @var nrhs         : number of right-hand-side terms
 * @var is_periodic  : periodic boundary condition is imposed or not
 * @var is_complex   : data type of the right-hand-side terms is fftw_complex or not (double)
 * @var l, c, u      : lower, center and upper-diagonal parts of the system
 * @var v            : internal buffer (updated "u" is stored)
 * @var q1           : internal buffer (additional right-hand-side term to be solved in addition to "q" is stored)
 */
struct tdm_info_t_ {
  int size;
  int nrhs;
  bool is_periodic;
  bool is_complex;
  double * restrict l;
  double * restrict c;
  double * restrict u;
  double * restrict v;
  double * restrict q1;
};

/**
 * @brief initialise tdm_info_t
 * @param[in]  size        : size of the system
 * @param[in]  nrhs        : number of right-hand-side terms
 * @param[in]  is_periodic : periodic boundary condition is imposed or not
 * @param[in]  is_complex  : data type of the right-hand-side terms is fftw_complex or not (double)
 * @param[out] info        : pointer to the resulting structure
 * @return                 : error code
 */
static int construct(
    const int size,
    const int nrhs,
    const bool is_periodic,
    const bool is_complex,
    tdm_info_t ** info
){
  // sanitise input
  if(size <= 0){
    printf("ERROR(%s): size should be positive: %d\n", __func__, size);
    *info = NULL;
    return 1;
  }
  if(nrhs <= 0){
    printf("ERROR(%s): nrhs should be positive: %d\n", __func__, nrhs);
    *info = NULL;
    return 1;
  }
  *info = memory_calloc(1, sizeof(tdm_info_t));
  // sizes
  (*info)->size = size;
  (*info)->nrhs = nrhs;
  // flags
  (*info)->is_periodic = is_periodic;
  (*info)->is_complex  = is_complex;
  // buffers
  (*info)->l = memory_calloc(size, sizeof(double));
  (*info)->c = memory_calloc(size, sizeof(double));
  (*info)->u = memory_calloc(size, sizeof(double));
  (*info)->v = memory_calloc(size, sizeof(double));
  if(is_periodic && /* to avoid zero-size allocation */ 1 < size){
    (*info)->q1 = memory_calloc(size - 1, sizeof(double));
  }else{
    (*info)->q1 = NULL;
  }
  return 0;
}

/**
 * @brief return pointer to the lower-diagonal matrix
 * @param[in]  info : initialised by constructor
 * @param[out] l    : pointer to the lower-diagonal matrix
 * @return          : error code
 */
static int get_l(
    const tdm_info_t * info,
    double * restrict * l
){
  if(NULL == info){
    printf("ERROR(%s): info is NULL\n", __func__);
    return 1;
  }
  *l = info->l;
  return 0;
}

/**
 * @brief return pointer to the center-diagonal matrix
 * @param[in]  info : initialised by constructor
 * @param[out] c    : pointer to the center-diagonal matrix
 * @return          : error code
 */
static int get_c(
    const tdm_info_t * info,
    double * restrict * c
){
  if(NULL == info){
    printf("ERROR(%s): info is NULL\n", __func__);
    return 1;
  }
  *c = info->c;
  return 0;
}

/**
 * @brief return pointer to the upper-diagonal matrix
 * @param[in]  info : initialised by constructor
 * @param[out] u    : pointer to the upper-diagonal matrix
 * @return          : error code
 */
static int get_u(
    const tdm_info_t * info,
    double * restrict * u
){
  if(NULL == info){
    printf("ERROR(%s): info is NULL\n", __func__);
    return 1;
  }
  *u = info->u;
  return 0;
}

static int get_size(
    const tdm_info_t * info,
    int * size
){
  if(NULL == info){
    printf("ERROR(%s): info is NULL\n", __func__);
    return 1;
  }
  *size = info->size;
  return 0;
}

static int get_nrhs(
    const tdm_info_t * info,
    int * nrhs
){
  if(NULL == info){
    printf("ERROR(%s): info is NULL\n", __func__);
    return 1;
  }
  *nrhs = info->nrhs;
  return 0;
}

/**
 * @brief solve tri-diagonal systems for the given input
 * @param[in]     info : initialised by constructor
 * @param[in,out] data : pointer to the right-hand-side terms, also used as a place to store the result
 * @return             : error code
 */
static int solve(
    const tdm_info_t * info,
    void * restrict data
){
  if(NULL == info){
    printf("ERROR(%s): info is NULL\n", __func__);
    return 1;
  }
  const int size = info->size;
  const int nrhs = info->nrhs;
  const bool is_periodic = info->is_periodic;
  const bool is_complex  = info->is_complex;
  const double * restrict l = info->l;
  const double * restrict c = info->c;
  const double * restrict u = info->u;
  double * restrict v  = info->v;
  double * restrict q1 = info->q1;
  if(is_complex){
    tdm_solve_fftw_complex(size, nrhs, is_periodic, l, c, u, v, data, q1);
  }else{
    tdm_solve_double(size, nrhs, is_periodic, l, c, u, v, data, q1);
  }
  return 0;
}

/**
 * @brief deallocate internal buffers
 * @param[in] info : initialised by constructor
 * @return         : error code
 */
static int destruct(
    tdm_info_t * info
){
  if(NULL == info){
    printf("ERROR(%s): info is NULL\n", __func__);
    return 1;
  }
  memory_free(info->l);
  memory_free(info->c);
  memory_free(info->u);
  memory_free(info->v);
  memory_free(info->q1);
  memory_free(info);
  return 0;
}

const tdm_t tdm = {
  .construct = construct,
  .get_l     = get_l,
  .get_c     = get_c,
  .get_u     = get_u,
  .get_size  = get_size,
  .get_nrhs  = get_nrhs,
  .solve     = solve,
  .destruct  = destruct,
};

