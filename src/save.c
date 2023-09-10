#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "sdecomp.h"
#include "param.h"
#include "memory.h"
#include "domain.h"
#include "save.h"
#include "fileio.h"
#include "config.h"

// parameters to specify directory name
static const char g_dirname_prefix[] = {"output/save/step"};
static const int g_dirname_ndigits = 10;

// name of directory
static char * g_dirname = NULL;
static size_t g_dirname_nchars = 0;

// scheduler
static double g_rate = 0.;
static double g_next = 0.;

/**
 * @brief constructor - schedule saving flow fields
 * @param[in] domain : MPI communicator
 * @param[in] time   : current time (hereafter in free-fall time units)
 */
static int init(
    const domain_t * domain,
    const double time
){
  // fetch timings
  if(0 != config.get_double("save_rate", &g_rate)){
    return 1;
  }
  double after = 0.;
  if(0 != config.get_double("save_after", &after)){
    return 1;
  }
  // schedule next event
  g_next = g_rate * ceil(
      fmax(DBL_EPSILON, fmax(time, after)) / g_rate
  );
  // allocate directory name
  g_dirname_nchars =
    + strlen(g_dirname_prefix)
    + g_dirname_ndigits;
  g_dirname = memory_calloc(g_dirname_nchars + 2, sizeof(char));
  // report
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  if(root == myrank){
    FILE * stream = stdout;
    fprintf(stream, "SAVE\n");
    fprintf(stream, "\tdest: %s\n", g_dirname_prefix);
    fprintf(stream, "\tnext: % .3e\n", g_next);
    fprintf(stream, "\trate: % .3e\n", g_rate);
    fflush(stream);
  }
  return 0;
}

/**
 * @brief prepare place to output flow fields
 * @param[in]  domain  : information related to MPI domain decomposition
 * @param[in]  step    : time step
 * @param[out] dirname : name of created directory
 */
static int prepare(
    const domain_t * domain,
    const int step,
    char ** dirname
){
  // set directory name
  snprintf(
      g_dirname,
      g_dirname_nchars + 1,
      "%s%0*d",
      g_dirname_prefix,
      g_dirname_ndigits,
      step
  );
  *dirname = g_dirname;
  // get communicator to identify the main process
  const int root = 0;
  int myrank = root;
  sdecomp.get_comm_rank(domain->info, &myrank);
  // create directory
  if(root == myrank){
    // although it may fail, anyway continue, which is designed to be safe
    fileio.mkdir(*dirname);
  }
  // wait for the main process to complete making directory
  MPI_Barrier(MPI_COMM_WORLD);
  // schedule next saving event
  g_next += g_rate;
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

const save_t save = {
  .init          = init,
  .prepare       = prepare,
  .get_next_time = get_next_time,
};

