#if !defined(CONFIG_H)
#define CONFIG_H

typedef struct {
  // getters for a double-precision value
  int (* const get_double)(
      const char dsetname[],
      double * value
  );
} config_t;

extern const config_t config;

#endif // CONFIG_H
