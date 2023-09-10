#include <stdio.h>
#include <math.h>
#include <float.h>
#include "memory.h"
#include "config.h"
#include "domain.h"
#include "fluid.h"
#include "fileio.h"
#include "logging.h"
#include "internal.h"

static double g_rate = 0.;
static double g_next = 0.;

/**
 * @brief constructor - schedule logging
 * @param[in] domain : MPI communicator
 * @param[in] time   : current time (hereafter in free-fall time units)
 */
static int init(
    const domain_t * domain,
    const double time
){
  if(0 != config.get_double("log_rate", &g_rate)){
    return 1;
  }
  g_next = g_rate * ceil(
      fmax(DBL_EPSILON, time) / g_rate
  );
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(root == myrank){
    printf("LOGGING\n");
    printf("\tnext: % .3e\n", g_next);
    printf("\trate: % .3e\n", g_rate);
    fflush(stdout);
  }
  return 0;
}

/**
 * @brief show current step, time, time step size, diffusive treatments
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] fname  : file name to which the log is written
 * @param[in] time   : current simulation time
 * @param[in] step   : current time step
 * @param[in] dt     : time step size
 * @param[in] wtime  : current wall time
 */
static void show_progress(
    const char fname[],
    const domain_t * domain,
    const double time,
    const size_t step,
    const double dt,
    const double wtime
){
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(root == myrank){
    FILE * fp = fileio.fopen(fname, "a");
    if(NULL != fp){
      // show progress to standard output and file | 8
      // output to stdout and file
#define MPRINT(...) { \
      fprintf(fp,     __VA_ARGS__); \
      fprintf(stdout, __VA_ARGS__); \
}
      MPRINT("step %zu, time %.1f, dt %.2e, elapsed %.1f [sec]\n", step, time, dt, wtime);
#undef MPRINT
      fileio.fclose(fp);
    }
  }
}

/**
 * @brief output log files to be monitored during simulation
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] step   : current time step
 * @param[in] time   : current simulation time
 * @param[in] dt     : time step size
 * @param[in] wtime  : current wall time
 * @param[in] fluid  : velocity and temperature
 */
static void check_and_output(
    const domain_t * domain,
    const size_t step,
    const double time,
    const double dt,
    const double wtime,
    const fluid_t * fluid
){
  show_progress           ("output/log/progress.dat",   domain, time, step, dt, wtime);
  logging_check_divergence("output/log/divergence.dat", domain, time, fluid);
  logging_check_momentum  ("output/log/momentum.dat",   domain, time, fluid);
  logging_check_energy    ("output/log/energy.dat",     domain, time, fluid);
  logging_check_nusselt   ("output/log/nusselt.dat",    domain, time, fluid);
  g_next += g_rate;
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

const logging_t logging = {
  .init             = init,
  .check_and_output = check_and_output,
  .get_next_time    = get_next_time,
};

