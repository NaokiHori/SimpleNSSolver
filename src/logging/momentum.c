#include <stdio.h>
#include "domain.h"
#include "fluid.h"
#include "fileio.h"
#include "array_macros/domain/dxf.h"
#include "array_macros/domain/dxc.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#if NDIMS == 3
#include "array_macros/fluid/uz.h"
#endif
#include "internal.h"

/**
 * @brief compute total momenta
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : information about domain decomposition and size
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity
 * @return           : error code
 */
int logging_check_momentum(
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
  const int isize = domain->mysizes[0];
  const int jsize = domain->mysizes[1];
#if NDIMS == 3
  const int ksize = domain->mysizes[2];
#endif
  const double * restrict dxf = domain->dxf;
  const double * restrict dxc = domain->dxc;
  const double dy = domain->dy;
#if NDIMS == 3
  const double dz = domain->dz;
#endif
  const double * restrict ux = fluid->ux.data;
  const double * restrict uy = fluid->uy.data;
#if NDIMS == 3
  const double * restrict uz = fluid->uz.data;
#endif
  double moms[NDIMS] = {0.};
  // compute total x-momentum | 19
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 2; i <= isize; i++){
      const double dx = DXC(i  );
      const double cellsize = dx * dy;
      moms[0] += UX(i, j) * cellsize;
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 2; i <= isize; i++){
        const double dx = DXC(i  );
        const double cellsize = dx * dy * dz;
        moms[0] += UX(i, j, k) * cellsize;
      }
    }
  }
#endif
  // compute total y-momentum | 19
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      const double dx = DXF(i  );
      const double cellsize = dx * dy;
      moms[1] += UY(i, j) * cellsize;
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        const double dx = DXF(i  );
        const double cellsize = dx * dy * dz;
        moms[1] += UY(i, j, k) * cellsize;
      }
    }
  }
#endif
#if NDIMS == 3
  // compute total z-momentum | 9
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        const double dx = DXF(i  );
        const double cellsize = dx * dy * dz;
        moms[2] += UZ(i, j, k) * cellsize;
      }
    }
  }
#endif
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : moms;
  void * recvbuf = moms;
  MPI_Reduce(sendbuf, recvbuf, NDIMS, MPI_DOUBLE, MPI_SUM, root, comm_cart);
  if(root == myrank){
    FILE * fp = fileio.fopen(fname, "a");
    if(NULL == fp){
      return 0;
    }
    fprintf(fp, "%8.2f ", time);
    for(int n = 0; n < NDIMS; n++){
      fprintf(fp, "% 18.15e%c", moms[n], NDIMS - 1 == n ? '\n' : ' ');
    }
    fileio.fclose(fp);
  }
  return 0;
}

