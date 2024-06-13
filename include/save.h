#if !defined(SAVE_H)
#define SAVE_H

#include "domain.h"
#include "fluid.h"

typedef struct save_t_ {
  // constructor
  int (* const init)(
    const domain_t * domain,
    const double time
  );
  // save a instantaneous flow field to files
  int (* const output)(
      const domain_t * domain,
      const size_t step,
      const double time,
      const fluid_t * fluid
  );
  // getter, next timing to call "output"
  double (* const get_next_time)(
      void
  );
} save_t;

extern const save_t save;

#endif // SAVE_H
