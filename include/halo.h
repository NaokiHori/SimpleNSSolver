#if !defined(HALO_H)
#define HALO_H

int halo_communicate_in_y(
    const domain_t * domain,
    MPI_Datatype * dtype,
    array_t * array
);

int halo_communicate_in_z(
    const domain_t * domain,
    MPI_Datatype * dtype,
    array_t * array
);

#endif // HALO_H
