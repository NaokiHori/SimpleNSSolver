#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "sdecomp.h"
#include "memory.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "statistics.h"
#include "fileio.h"
#include "config.h"
#include "array_macros/domain/hxxf.h"
#include "array_macros/domain/jdxf.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#include "array_macros/fluid/uz.h"
#include "array_macros/fluid/t.h"
#include "array_macros/statistics/ux1.h"
#include "array_macros/statistics/ux2.h"
#include "array_macros/statistics/uy1.h"
#include "array_macros/statistics/uy2.h"
#include "array_macros/statistics/uz1.h"
#include "array_macros/statistics/uz2.h"
#include "array_macros/statistics/t1.h"
#include "array_macros/statistics/t2.h"
#include "array_macros/statistics/adv.h"
#include "array_macros/statistics/dif.h"

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
static array_t g_uz1 = {0};
static array_t g_uz2 = {0};
static array_t g_t1 = {0};
static array_t g_t2 = {0};
static array_t g_adv = {0};
static array_t g_dif = {0};

/**
 * @brief constructor - initialise and allocate internal buffers, schedule collection
 * @param[in] domain : information about domain decomposition and size
 * @param[in] time   : current time (hereafter in free-fall time units)
 * @return           : error code
 */
static int init (
    const domain_t * domain,
    const double time
) {
  // fetch timings
  if (0 != config.get_double("stat_rate", &g_rate)) {
    return 1;
  }
  double after = 0.;
  if (0 != config.get_double("stat_after", &after)) {
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
  if (0 != array.create(domain, UX1_NADDS, sizeof(double), &g_ux1)) return 1;
  if (0 != array.create(domain, UX2_NADDS, sizeof(double), &g_ux2)) return 1;
  if (0 != array.create(domain, UY1_NADDS, sizeof(double), &g_uy1)) return 1;
  if (0 != array.create(domain, UY2_NADDS, sizeof(double), &g_uy2)) return 1;
  if (0 != array.create(domain, UZ1_NADDS, sizeof(double), &g_uz1)) return 1;
  if (0 != array.create(domain, UZ2_NADDS, sizeof(double), &g_uz2)) return 1;
  if (0 != array.create(domain, T1_NADDS,  sizeof(double), &g_t1 )) return 1;
  if (0 != array.create(domain, T2_NADDS,  sizeof(double), &g_t2 )) return 1;
  if (0 != array.create(domain, ADV_NADDS, sizeof(double), &g_adv)) return 1;
  if (0 != array.create(domain, DIF_NADDS, sizeof(double), &g_dif)) return 1;
  // report
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if (root == myrank) {
    FILE * stream = stdout;
    fprintf(stream, "STATISTICS\n");
    fprintf(stream, "\tdest: %s\n", g_dirname_prefix);
    fprintf(stream, "\tnext: % .3e\n", g_next);
    fprintf(stream, "\trate: % .3e\n", g_rate);
    fflush(stream);
  }
  return 0;
}

static double get_next_time (
    void
) {
  return g_next;
}

static void collect_mean_ux (
    const domain_t * domain,
    const double * restrict ux
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  double * restrict ux1 = g_ux1.data;
  double * restrict ux2 = g_ux2.data;
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize + 1; i++) {
        UX1(i, j, k) += pow(UX(i, j, k), 1.);
        UX2(i, j, k) += pow(UX(i, j, k), 2.);
      }
    }
  }
}

static void collect_mean_uy (
    const domain_t * domain,
    const double * restrict uy
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  double * restrict uy1 = g_uy1.data;
  double * restrict uy2 = g_uy2.data;
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 0; i <= isize + 1; i++) {
        UY1(i, j, k) += pow(UY(i, j, k), 1.);
        UY2(i, j, k) += pow(UY(i, j, k), 2.);
      }
    }
  }
}

static void collect_mean_uz (
    const domain_t * domain,
    const double * restrict uz
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  double * restrict uz1 = g_uz1.data;
  double * restrict uz2 = g_uz2.data;
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 0; i <= isize + 1; i++) {
        UZ1(i, j, k) += pow(UZ(i, j, k), 1.);
        UZ2(i, j, k) += pow(UZ(i, j, k), 2.);
      }
    }
  }
}

static void collect_mean_t (
    const domain_t * domain,
    const double * restrict t
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  double * restrict t1 = g_t1.data;
  double * restrict t2 = g_t2.data;
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 0; i <= isize + 1; i++) {
        T1(i, j, k) += pow(T(i, j, k), 1.);
        T2(i, j, k) += pow(T(i, j, k), 2.);
      }
    }
  }
}

static void collect_adv (
    const domain_t * domain,
    const double * restrict ux,
    const double * restrict t
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  double * restrict adv = g_adv.data;
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize + 1; i++) {
        ADV(i, j, k) += JDXF(i) / HXXF(i) * UX(i, j, k) * (0.5 * T(i-1, j  , k  ) + 0.5 * T(i  , j  , k  ));
      }
    }
  }
}

static void collect_dif (
    const domain_t * domain,
    const double diffusivity,
    const double * restrict t
) {
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
  const int ksize = domain->mysizes[2];
  const double * restrict hxxf = domain->hxxf;
  const double * restrict jdxf = domain->jdxf;
  double * restrict dif = g_dif.data;
  for (int k = 1; k <= ksize; k++) {
    for (int j = 1; j <= jsize; j++) {
      for (int i = 1; i <= isize + 1; i++) {
        const double t_xm = T(i-1, j  , k  );
        const double t_xp = T(i  , j  , k  );
        DIF(i, j, k) -= diffusivity * JDXF(i) / HXXF(i) / HXXF(i) * (t_xp - t_xm);
      }
    }
  }
}

/**
 * @brief accumulate statistical data
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] fluid  : flow field
 * @return           : error code
 */
static int collect (
    const domain_t * domain,
    const fluid_t * fluid
) {
  // collect temporally-averaged quantities
  collect_mean_ux(domain, fluid->ux.data);
  collect_mean_uy(domain, fluid->uy.data);
  collect_mean_uz(domain, fluid->uz.data);
  collect_mean_t(domain, fluid->t.data);
  collect_adv(domain, fluid->ux.data, fluid->t.data);
  collect_dif(domain, fluid_compute_temperature_diffusivity(fluid), fluid->t.data);
  // increment number of samples
  g_num += 1;
  // schedule next event
  g_next += g_rate;
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
) {
  // when no statistics are collected (g_num is 0),
  //   no reason to save, so abort
  if (0 == g_num) {
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
  if (root == myrank) {
    // although it may fail, anyway continue, which is designed to be safe
    fileio.mkdir(g_dirname);
    // save scalars
    fileio.w_serial(g_dirname, "num", 0, NULL, fileio.npy_size_t, sizeof(size_t), &g_num);
  }
  // wait for the main process to complete making directory
  MPI_Barrier(MPI_COMM_WORLD);
  // save domain info (coordinates)
  domain_save(g_dirname, domain);
  // save collected statistics
  // prepare list of all variables to be dumped
  typedef struct {
    const char * name;
    const array_t * array;
  } variable_t;
  const variable_t variables[] = {
    {.name = "ux1", .array = &g_ux1},
    {.name = "ux2", .array = &g_ux2},
    {.name = "uy1", .array = &g_uy1},
    {.name = "uy2", .array = &g_uy2},
    {.name = "uz1", .array = &g_uz1},
    {.name = "uz2", .array = &g_uz2},
    {.name = "t1",  .array = &g_t1 },
    {.name = "t2",  .array = &g_t2 },
    {.name = "adv", .array = &g_adv},
    {.name = "dif", .array = &g_dif},
  };
  // dump each statistical variable to a file
  for (size_t n = 0; n < sizeof(variables) / sizeof(variables[0]); n++) {
    const variable_t * v = variables + n;
    array.dump(domain, g_dirname, v->name, fileio.npy_double, v->array);
  }
  return 0;
}

const statistics_t statistics = {
  .init          = init,
  .collect       = collect,
  .output        = output,
  .get_next_time = get_next_time,
};

