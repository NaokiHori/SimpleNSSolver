#include <stdio.h>
#include <mpi.h>
#include "memory.h"
#include "timer.h"
#include "domain.h"
#include "fluid.h"
#include "fluid_solver.h"
#include "integrate.h"
#include "statistics.h"
#include "save.h"
#include "logging.h"
#include "config.h"
#include "fileio.h"

/**
 * @brief main function
 * @param[in] argc : number of arguments (expect 2)
 * @param[in] argv : name of the directory
 *                     where a set of initial condition is contained
 * @return         : error code
 */
int main (
    int argc,
    char * argv[]
) {
  // launch MPI, start timer | 5
  MPI_Init(NULL, NULL);
  const int root = 0;
  int myrank = root;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  const double tic = timer();
  // find name of directory where IC is stored | 7
  if (2 != argc) {
    if (root == myrank) {
      printf("directory name should be given as input\n");
    }
    goto abort;
  }
  const char * dirname_ic = argv[1];
  // initialise fileio object
  if (0 != fileio.init()) {
    goto abort;
  }
  // initialise time step and time units | 8
  size_t step = 0;
  if (0 != fileio.r_serial(dirname_ic, "step", 0, NULL, fileio.npy_size_t, sizeof(size_t), &step)) {
    goto abort;
  }
  double time = 0.;
  if (0 != fileio.r_serial(dirname_ic, "time", 0, NULL, fileio.npy_double, sizeof(double), &time)) {
    goto abort;
  }
  // initialise structures | 8
  domain_t domain = {0};
  if (0 != domain_init(dirname_ic, &domain)) {
    goto abort;
  }
  fluid_t fluid = {0};
  if (0 != fluid_init(dirname_ic, &domain, &fluid)) {
    goto abort;
  }
  // initialise auxiliary objects | 9
  if (0 != logging.init(&domain, time)) {
    goto abort;
  }
  if (0 != save.init(&domain, time)) {
    goto abort;
  }
  if (0 != statistics.init(&domain, time)) {
    goto abort;
  }
  // check termination conditions
  double timemax = 0.;
  if (0 != config.get_double("timemax", &timemax)) {
    goto abort;
  }
  double wtimemax = 0.;
  if (0 != config.get_double("wtimemax", &wtimemax)) {
    goto abort;
  }
  // report
  if (root == myrank) {
    printf("step: %zu, time: % .7e\n", step, time);
    printf("timemax: % .7e, wtimemax: % .7e\n", timemax, wtimemax);
  }
  // main loop
  for (;;) {
    // proceed for one step
    double dt = 0.;
    if (0 != integrate(&domain, &fluid, &dt)) {
      goto abort;
    }
    // update step and simulation / wall time | 3
    step += 1;
    time += dt;
    const double toc = timer();
    // terminate if one of the following conditions is met | 8
    // the simulation is finished
    if (timemax < time) {
      break;
    }
    // wall time limit is reached
    if (wtimemax < toc - tic) {
      break;
    }
    // compute and output log regularly | 3
    if (logging.get_next_time() < time) {
      logging.check_and_output(&domain, step, time, dt, toc - tic, &fluid);
    }
    // save flow fields regularly | 3
    if (save.get_next_time() < time) {
      save.output(&domain, step, time, &fluid);
    }
    // collect statistics regularly | 3
    if (statistics.get_next_time() < time) {
      statistics.collect(&domain, &fluid);
    }
  }
  // finalisation
  // save final flow fields | 1
  save.output(&domain, step, time, &fluid);
  // save collected statistics | 1
  statistics.output(&domain, step);
  // finalise MPI | 2
abort:
  MPI_Finalize();
  return 0;
}

