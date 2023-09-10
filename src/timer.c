#include <mpi.h>
#include "timer.h"

/**
 * @brief get current time
 * @return : current time
 */
double timer(
    void
){
  const int root = 0;
  // although this is called by all processes,
  double time = MPI_Wtime();
  // use the result of the main process
  MPI_Bcast(&time, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  return time;
}

