#include "param.h"

const bool param_m_implicit_x =  true;
const bool param_m_implicit_y = false;
#if NDIMS == 3
const bool param_m_implicit_z = false;
#endif

const bool param_t_implicit_x =  true;
const bool param_t_implicit_y = false;
#if NDIMS == 3
const bool param_t_implicit_z = false;
#endif

