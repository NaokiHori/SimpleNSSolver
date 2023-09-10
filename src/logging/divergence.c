#include <stdio.h>
#include <math.h>
#include "domain.h"
#include "fluid.h"
#include "fileio.h"
#include "array_macros/domain/dxf.h"
#include "array_macros/fluid/ux.h"
#include "array_macros/fluid/uy.h"
#if NDIMS == 3
#include "array_macros/fluid/uz.h"
#endif
#include "internal.h"

/**
 * @brief check divergence and write the maximum value
 * @param[in] fname  : file name to which the log is written
 * @param[in] domain : domain information
 * @param[in] time   : current simulation time
 * @param[in] fluid  : velocity
 * @return           : error code
 */
int logging_check_divergence(
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
  const double dy = domain->dy;
#if NDIMS == 3
  const double dz = domain->dz;
#endif
  const double * restrict ux = fluid->ux.data;
  const double * restrict uy = fluid->uy.data;
#if NDIMS == 3
  const double * restrict uz = fluid->uz.data;
#endif
  double divmax = 0.;
#if NDIMS == 2
  for(int j = 1; j <= jsize; j++){
    for(int i = 1; i <= isize; i++){
      // compute local divergence | 8
      const double dx = DXF(i  );
      const double ux_xm = UX(i  , j  );
      const double ux_xp = UX(i+1, j  );
      const double uy_ym = UY(i  , j  );
      const double uy_yp = UY(i  , j+1);
      const double div =
        +(ux_xp - ux_xm) / dx
        +(uy_yp - uy_ym) / dy;
      // check maximum | 1
      divmax = fmax(divmax, fabs(div));
    }
  }
#else
  for(int k = 1; k <= ksize; k++){
    for(int j = 1; j <= jsize; j++){
      for(int i = 1; i <= isize; i++){
        // compute local divergence | 11
        const double dx = DXF(i  );
        const double ux_xm = UX(i  , j  , k  );
        const double ux_xp = UX(i+1, j  , k  );
        const double uy_ym = UY(i  , j  , k  );
        const double uy_yp = UY(i  , j+1, k  );
        const double uz_zm = UZ(i  , j  , k  );
        const double uz_zp = UZ(i  , j  , k+1);
        const double div =
          +(ux_xp - ux_xm) / dx
          +(uy_yp - uy_ym) / dy
          +(uz_zp - uz_zm) / dz;
        // check maximum | 1
        divmax = fmax(divmax, fabs(div));
      }
    }
  }
#endif
  // collect information among all processes
  const void * sendbuf = root == myrank ? MPI_IN_PLACE : &divmax;
  void * recvbuf = &divmax;
  MPI_Reduce(sendbuf, recvbuf, 1, MPI_DOUBLE, MPI_MAX, root, comm_cart);
  // result is written to a file from the main process
  if(root == myrank){
    FILE * fp = fileio.fopen(fname, "a");
    if(NULL == fp){
      return 0;
    }
    fprintf(fp, "%8.2f % .1e\n", time, divmax);
    fileio.fclose(fp);
  }
  return 0;
}

