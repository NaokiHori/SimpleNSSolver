#if !defined(TDM_H)
#define TDM_H

#include <stdbool.h>

typedef struct tdm_info_t_ tdm_info_t;

typedef struct {
  int (* const construct)(
      const int size,
      const int nrhs,
      const bool is_periodic,
      const bool is_complex,
      tdm_info_t ** info
  );
  int (* const get_l)(
      const tdm_info_t * info,
      double * restrict * l
  );
  int (* const get_c)(
      const tdm_info_t * info,
      double * restrict * c
  );
  int (* const get_u)(
      const tdm_info_t * info,
      double * restrict * u
  );
  int (* const get_size)(
      const tdm_info_t * info,
      int * size
  );
  int (* const get_nrhs)(
      const tdm_info_t * info,
      int * nrhs
  );
  int (* const solve)(
      const tdm_info_t * info,
      void * restrict data
  );
  int (* const destruct)(
      tdm_info_t * info
  );
} tdm_t;

extern const tdm_t tdm;

#endif // TDM_H
