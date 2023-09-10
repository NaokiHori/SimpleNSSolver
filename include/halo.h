#if !defined(HALO_H)
#define HALO_H

int halo_communicate_in_y(
    const domain_t * domain,
    MPI_Datatype * dtype,
    array_t * array
);

#if NDIMS == 3
int halo_communicate_in_z(
    const domain_t * domain,
    MPI_Datatype * dtype,
    array_t * array
);
#endif

#endif // HALO_H
