#if !defined(INCLUDE_ARRAY_MACROS_STATISTICS_ADV_H)
#define INCLUDE_ARRAY_MACROS_STATISTICS_ADV_H

// This file is generated by tools/define_arrays.py

// [1 : isize+1], [1 : jsize+0]
#define ADV(I, J) (adv[(I-1) + (isize+1) * (J-1)])
#define ADV_NADDS (int [NDIMS][2]){ {0, 1}, {0, 0}, }

#endif // INCLUDE_ARRAY_MACROS_STATISTICS_ADV_H
