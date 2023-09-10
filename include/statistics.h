#if !defined(STATISTICS_H)
#define STATISTICS_H

#include "domain.h"
#include "fluid.h"

typedef struct {
  // constructor
  int (* const init)(
      const domain_t * domain,
      const double time
  );
  // collecting statistics
  int (* const collect)(
      const domain_t * domain,
      const fluid_t * fluid
  );
  // save statistics to files
  int (* const output)(
      const domain_t * domain,
      const size_t step
  );
  // getter, next timing to call "collect"
  double (* const get_next_time)(
      void
  );
} statistics_t;

extern const statistics_t statistics;

#endif // STATISTICS_H
