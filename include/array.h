#if !defined(ARRAY_H)
#define ARRAY_H

// struct and methods for multi-dimensional arrays
//   distributed among multiple processes

#include <stddef.h>
#include "domain.h"

typedef struct {
  // size of each element
  size_t size;
  // number of additional cells w.r.t. no-halo array
  //   i.e. 0th elements: [1 : mysizes[0]]
  //        1st elements: [1 : mysizes[1]]
  //        2nd elements: [1 : mysizes[2]]
  //        ...
  int (* nadds)[2];
  // total size of local array (i.e. product of mysizes)
  size_t datasize;
  // pointer to the raw local array
  void * data;
} array_t;

typedef struct {
  // allocate array and store its size information
  int (* const prepare)(
      const domain_t * domain,
      const int nadds[NDIMS][2],
      const size_t size,
      array_t * array
  );
  // clean-up local memory to store the array
  int (* const destroy)(
      array_t * array
  );
  // load array from NPY file
  int (* const load)(
      const domain_t * domain,
      const char dirname[],
      const char dsetname[],
      const char dtype[],
      array_t * array
  );
  // save array to NPY file
  int (* const dump)(
      const domain_t * domain,
      const char dirname[],
      const char dsetname[],
      const char dtype[],
      const array_t * array
  );
} array_method_t;

extern const array_method_t array;

#endif // ARRAY_H
