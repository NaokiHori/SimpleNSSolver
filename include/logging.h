#if !defined(LOGGING_H)
#define LOGGING_H

#include "domain.h"
#include "fluid.h"

typedef struct {
  // constructor
  int (* const init)(
      const domain_t * domain,
      const double time
  );
  // check quantities and dump to log files
  void (* const check_and_output)(
      const domain_t * domain,
      const size_t step,
      const double time,
      const double dt,
      const double wtime,
      const fluid_t * fluid
  );
  // getter, next timing to call "check_and_output"
  double (* const get_next_time)(
      void
  );
} logging_t;

extern const logging_t logging;

#endif // LOGGING_H
