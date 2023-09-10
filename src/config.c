#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <mpi.h>
#include "config.h"

/**
 * @brief load environment variable and interpret it as an double-precision value
 * @param[in]  dsetname : name of the environment variable
 * @param[out] value    : resulting value
 * @return              : error code
 */
static int get_double(
    const char dsetname[],
    double * value
){
  // error code
  int retval = 0;
  // try to load, may fail if the variable is not defined
  char * string = getenv(dsetname);
  if(NULL == string){
    retval = 1;
  }
  MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(0 != retval){
    printf("%s not found\n", dsetname);
    return 1;
  }
  // try to convert a string to a double-precision value
  errno = 0;
  *value = strtod(string, NULL);
  if(0 != errno){
    retval = 1;
  }
  MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(0 != retval){
    printf("%s: invalid value as double\n", dsetname);
    return 1;
  }
  return 0;
}

const config_t config = {
  .get_double = get_double,
};

