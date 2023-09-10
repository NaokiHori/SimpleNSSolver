#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "sdecomp.h"
#include "param.h"
#include "memory.h"
#include "domain.h"
#include "statistics.h"
#include "fileio.h"
#include "config.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#if NDIMS == 3
#include "array_macros/fluid/uz.h"
#endif
#include "array_macros/fluid/t.h"
#include "array_macros/statistics/ux1.h"
#include "array_macros/statistics/ux2.h"
#include "array_macros/statistics/uy1.h"
#include "array_macros/statistics/uy2.h"
#if NDIMS == 3
#include "array_macros/statistics/uz1.h"
#include "array_macros/statistics/uz2.h"
#endif
#include "array_macros/statistics/t1.h"
#include "array_macros/statistics/t2.h"
#include "array_macros/statistics/uxt.h"

// compress 3D statistical field to 1D (true)
//   or save it as it is (false)
static const bool g_reduction = true;

// parameters to specify directory name
static const char g_dirname_prefix[] = {"output/stat/step"};
static const int g_dirname_ndigits = 10;

// name of directory
static char * g_dirname = NULL;
static size_t g_dirname_nchars = 0;

// scheduler
static double g_rate = 0.;
static double g_next = 0.;

// data
static size_t g_num = 0;
static array_t g_ux1 = {0};
static array_t g_ux2 = {0};
static array_t g_uy1 = {0};
static array_t g_uy2 = {0};
#if NDIMS == 3
static array_t g_uz1 = {0};
static array_t g_uz2 = {0};
#endif
static array_t g_t1 = {0};
static array_t g_t2 = {0};
static array_t g_uxt = {0};

// diffusivities, store here for convenience
static double g_m_dif = 0.;
static double g_t_dif = 0.;

/**
 * @brief constructor - initialise and allocate internal buffers, schedule collection
 * @param[in] domain : information about domain decomposition and size
 * @param[in] time   : current time (hereafter in free-fall time units)
 * @return           : error code
 */
static int init(
    const domain_t * domain,
    const double time
){
  // fetch timings
  if(0 != config.get_double("stat_rate", &g_rate)){
    return 1;
  }
  double after = 0.;
  if(0 != config.get_double("stat_after", &after)){
    return 1;
  }
  g_next = g_rate * ceil(
      fmax(DBL_EPSILON, fmax(time, after)) / g_rate
  );
  // allocate directory name
  g_dirname_nchars =
    + strlen(g_dirname_prefix)
    + g_dirname_ndigits;
  g_dirname = memory_calloc(g_dirname_nchars + 2, sizeof(char));
  // prepare arrays
  if(0 != array.prepare(domain, UX1_NADDS, sizeof(double), &g_ux1)) return 1;
  if(0 != array.prepare(domain, UX2_NADDS, sizeof(double), &g_ux2)) return 1;
  if(0 != array.prepare(domain, UY1_NADDS, sizeof(double), &g_uy1)) return 1;
  if(0 != array.prepare(domain, UY2_NADDS, sizeof(double), &g_uy2)) return 1;
#if NDIMS == 3
  if(0 != array.prepare(domain, UZ1_NADDS, sizeof(double), &g_uz1)) return 1;
  if(0 != array.prepare(domain, UZ2_NADDS, sizeof(double), &g_uz2)) return 1;
#endif
  if(0 != array.prepare(domain, T1_NADDS,  sizeof(double), &g_t1 )) return 1;
  if(0 != array.prepare(domain, T2_NADDS,  sizeof(double), &g_t2 )) return 1;
  if(0 != array.prepare(domain, UXT_NADDS, sizeof(double), &g_uxt)) return 1;
  // report
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(root == myrank){
    FILE * stream = stdout;
    fprintf(stream, "STATISTICS\n");
    fprintf(stream, "\tdest: %s\n", g_dirname_prefix);
    fprintf(stream, "\tnext: % .3e\n", g_next);
    fprintf(stream, "\trate: % .3e\n", g_rate);
    fflush(stream);
  }
  return 0;
}

/**
 * @brief getter of a member: g_next
 * @return : g_next
 */
static double get_next_time(
    void
){
  return g_next;
}

/**
 * @brief compute ux^1 and ux^2 and add results to the arrays
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] ux     : x velocity
 */
static void collect_mean_ux(
    const domain_t * domain,
    const double * restrict ux
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  double * restrict ux1 = g_ux1.data;
  double * restrict ux2 = g_ux2.data;
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize + 1; i++){
      UX1(i, j) += pow(UX(i, j), 1.);
      UX2(i, j) += pow(UX(i, j), 2.);
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize + 1; i++){
        UX1(i, j, k) += pow(UX(i, j, k), 1.);
        UX2(i, j, k) += pow(UX(i, j, k), 2.);
      }
    }
  }
#endif
}

/**
 * @brief compute uy^1 and uy^2 and add results to the arrays
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] uy     : y velocity
 */
static void collect_mean_uy(
    const domain_t * domain,
    const double * restrict uy
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  double * restrict uy1 = g_uy1.data;
  double * restrict uy2 = g_uy2.data;
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 0; i <= isize + 1; i++){
      UY1(i, j) += pow(UY(i, j), 1.);
      UY2(i, j) += pow(UY(i, j), 2.);
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 0; i <= isize + 1; i++){
        UY1(i, j, k) += pow(UY(i, j, k), 1.);
        UY2(i, j, k) += pow(UY(i, j, k), 2.);
      }
    }
  }
#endif
}

#if NDIMS == 3
/**
 * @brief compute uz^1 and uz^2 and add results to the arrays
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] uz     : z velocity
 */
static void collect_mean_uz(
    const domain_t * domain,
    const double * restrict uz
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  double * restrict uz1 = g_uz1.data;
  double * restrict uz2 = g_uz2.data;
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 0; i <= isize + 1; i++){
        UZ1(i, j, k) += pow(UZ(i, j, k), 1.);
        UZ2(i, j, k) += pow(UZ(i, j, k), 2.);
      }
    }
  }
}
#endif

/**
 * @brief compute T^1 and T^2 and add results to the arrays
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] t      : temperature
 */
static void collect_mean_t(
    const domain_t * domain,
    const double * restrict t
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  double * restrict t1 = g_t1.data;
  double * restrict t2 = g_t2.data;
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 0; i <= isize + 1; i++){
      T1(i, j) += pow(T(i, j), 1.);
      T2(i, j) += pow(T(i, j), 2.);
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 0; i <= isize + 1; i++){
        T1(i, j, k) += pow(T(i, j, k), 1.);
        T2(i, j, k) += pow(T(i, j, k), 2.);
      }
    }
  }
#endif
}

/**
 * @brief compute ux T and add results to the arrays
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] ux     : x velocity
 * @param[in] t      : temperature
 */
static void collect_uxt(
    const domain_t * domain,
    const double * restrict ux,
    const double * restrict t
){
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  double * restrict uxt = g_uxt.data;
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize + 1; i++){
      const double t_ =
        + 0.5 * T(i-1, j  )
        + 0.5 * T(i  , j  );
      UXT(i, j) += UX(i, j) * t_;
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize + 1; i++){
        const double t_ =
          + 0.5 * T(i-1, j  , k  )
          + 0.5 * T(i  , j  , k  );
        UXT(i, j, k) += UX(i, j, k) * t_;
      }
    }
  }
#endif
}

/**
 * @brief accumulate statistical data
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] fluid  : flow field
 * @return           : error code
 */
static int collect(
    const domain_t * domain,
    const fluid_t * fluid
){
  // collect temporally-averaged quantities
  collect_mean_ux(domain, fluid->ux.data);
  collect_mean_uy(domain, fluid->uy.data);
#if NDIMS == 3
  collect_mean_uz(domain, fluid->uz.data);
#endif
  collect_mean_t(domain, fluid->t.data);
  collect_uxt(domain, fluid->ux.data, fluid->t.data);
  // assign diffusivities
  g_m_dif = fluid->m_dif;
  g_t_dif = fluid->t_dif;
  // increment number of samples
  g_num += 1;
  // schedule next event
  g_next += g_rate;
  return 0;
}

/**
 * @brief reduce N-dimensional (xy, xyz) array to 1D (x) vector
 * @param[in] domain   : information related to MPI domain decomposition
 * @param[in] dirname  : name of directory
 * @param[in] dsetname : name of dataset
 * @param[in] array    : N-dimensional array to be reduced
 */
static int reduce_and_write(
    const domain_t * domain,
    const char dirname[],
    const char dsetname[],
    const array_t * array
){
  // NOTE: assuming no halo cells in y
  if(0 != array->nadds[1][0] || 0 != array->nadds[1][1]){
    printf("stat array seems to have irregular shape in y\n");
    return 1;
  }
#if NDIMS == 3
  // NOTE: assuming no halo cells in z
  if(0 != array->nadds[2][0] || 0 != array->nadds[2][1]){
    printf("stat array seems to have irregular shape in z\n");
    return 1;
  }
#endif
  const int isize = domain->mysizes[0]
                  + array->nadds[0][0]
                  + array->nadds[0][1];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict data = array->data;
  double * restrict vec = memory_calloc(isize, sizeof(double));
#if NDIMS == 2
  for(int j = 0; j < jsize; j++){
    for(int i = 0; i < isize; i++){
      vec[i] += data[j * isize + i];
    }
  }
#else
  for(int k = 0; k < ksize; k++){
    for(int j = 0; j < jsize; j++){
      for(int i = 0; i < isize; i++){
        vec[i] += data[k * jsize * isize + j * isize + i];
      }
    }
  }
#endif
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : vec;
  void * recvbuf = vec;
  MPI_Reduce(sendbuf, recvbuf, isize, MPI_DOUBLE, MPI_SUM, root, comm_cart);
  if(root == myrank){
    fileio.w_serial(
        dirname,
        dsetname,
        1,
        (size_t [1]){isize},
        fileio.npy_double,
        sizeof(double),
        vec
    );
  }
  memory_free(vec);
  return 0;
}

/**
 * @brief save structures which contains collected statistical data
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] step   : current time step
 * @return           : error code
 */
static int output(
    const domain_t * domain,
    const size_t step
){
  // when no statistics are collected (g_num is 0),
  //   no reason to save, so abort
  if(0 == g_num){
    return 0;
  }
  // set directory name
  snprintf(
      g_dirname,
      g_dirname_nchars + 1,
      "%s%0*zu",
      g_dirname_prefix,
      g_dirname_ndigits,
      step
  );
  // get communicator to identify the main process
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  // create directory and save scalars from main process
  if(root == myrank){
    // although it may fail, anyway continue, which is designed to be safe
    fileio.mkdir(g_dirname);
    // save scalars
    fileio.w_serial(g_dirname, "num", 0, NULL, fileio.npy_size_t, sizeof(size_t), &g_num);
    fileio.w_serial(g_dirname, "m_dif", 0, NULL, fileio.npy_double, sizeof(double), &g_m_dif);
    fileio.w_serial(g_dirname, "t_dif", 0, NULL, fileio.npy_double, sizeof(double), &g_t_dif);
  }
  // wait for the main process to complete making directory
  MPI_Barrier(MPI_COMM_WORLD);
  // save domain info (coordinates)
  domain_save(g_dirname, domain);
  // save collected statistics
  if(g_reduction){
    reduce_and_write(domain, g_dirname, "ux1", &g_ux1);
    reduce_and_write(domain, g_dirname, "ux2", &g_ux2);
    reduce_and_write(domain, g_dirname, "uy1", &g_uy1);
    reduce_and_write(domain, g_dirname, "uy2", &g_uy2);
#if NDIMS == 3
    reduce_and_write(domain, g_dirname, "uz1", &g_uz1);
    reduce_and_write(domain, g_dirname, "uz2", &g_uz2);
#endif
    reduce_and_write(domain, g_dirname,  "t1", &g_t1);
    reduce_and_write(domain, g_dirname,  "t2", &g_t2);
    reduce_and_write(domain, g_dirname, "uxt", &g_uxt);
  }else{
    array.dump(domain, g_dirname, "ux1", fileio.npy_double, &g_ux1);
    array.dump(domain, g_dirname, "ux2", fileio.npy_double, &g_ux2);
    array.dump(domain, g_dirname, "uy1", fileio.npy_double, &g_uy1);
    array.dump(domain, g_dirname, "uy2", fileio.npy_double, &g_uy2);
#if NDIMS == 3
    array.dump(domain, g_dirname, "uz1", fileio.npy_double, &g_uz1);
    array.dump(domain, g_dirname, "uz2", fileio.npy_double, &g_uz2);
#endif
    array.dump(domain, g_dirname,  "t1", fileio.npy_double, &g_t1);
    array.dump(domain, g_dirname,  "t2", fileio.npy_double, &g_t2);
    array.dump(domain, g_dirname, "uxt", fileio.npy_double, &g_uxt);
  }
  return 0;
}

const statistics_t statistics = {
  .init          = init,
  .collect       = collect,
  .output        = output,
  .get_next_time = get_next_time,
};

