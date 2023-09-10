#include <stdio.h>
#include <mpi.h>
#include "domain.h"
#include "fluid.h"
#include "fileio.h"
#include "internal.h"

/**
 * @brief compute Nusselt number based on various definitions
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : information related to MPI domain decomposition
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity, temperature, and their diffusivities
 * @return           : error code
 */
int logging_check_nusselt(
    const char fname[],
    const domain_t * domain,
    const double time,
    const fluid_t * fluid
){
  const int root = 0;
  int myrank = root;
  MPI_Comm comm_cart = MPI_COMM_NULL;
  sdecomp.get_comm_rank(domain->info, &myrank);
  sdecomp.get_comm_cart(domain->info, &comm_cart);
  // compute Nusselt in several ways
  const double results[] = {
    // compute heat flux on the walls
    logging_internal_compute_nu_heat_flux(domain, fluid),
    // compute kinetic energy injection
    logging_internal_compute_nu_kinetic_energy_injection(domain, fluid),
    // compute kinetic energy dissipation
    logging_internal_compute_nu_kinetic_energy_dissipation(domain, fluid),
    // comoute thermal energy dissipation
    logging_internal_compute_nu_thermal_energy_dissipation(domain, fluid),
  };
  if(root == myrank){
    FILE * fp = fileio.fopen(fname, "a");
    if(NULL == fp){
      return 0;
    }
    fprintf(fp, "%8.2f ", time);
    const int ntypes = sizeof(results) / sizeof(results[0]);
    for(int n = 0; n < ntypes; n++){
      fprintf(fp, "% 18.15e%c", results[n], ntypes - 1 == n ? '\n' : ' ');
    }
    fileio.fclose(fp);
  }
  return 0;
}

