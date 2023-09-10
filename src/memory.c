#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "memory.h"

/**
 * @brief general-purpose memory allocator
 * @param[in] count : number of elements to be allocated
 * @param[in] size  : size of each element
 * @return          : pointer to the allocated buffer
 */
void * memory_calloc(
    const size_t count,
    const size_t size
){
  // try to allocate
  void * ptr = calloc(count, size);
  if(NULL == ptr){
    // failed to allocate
    FILE * stream = stderr;
    fprintf(stream, "memory allocation error: calloc(%zu, %zu)\n", count, size);
    fflush(stream);
    // since memory errors are fatal, I abort the program
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  return ptr;
}

/**
 * @brief corresponding memory deallocator
 * @param[in] ptr : pointer to the allocated buffer
 */
void memory_free(
    void * ptr
){
  // for now just wrap free
  free(ptr);
}

