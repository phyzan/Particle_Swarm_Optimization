#ifndef __FUNC__
#define __FUNC__

#include "../local_definitions.hpp"

double uniform_real(double a, double b);
double gauss_distribution(double mean, double stddev);

MPFR_ARR double_to_mpfr_ARR(DOUBLE_ARR arr);
MPFR_VEC double_to_mpfr_VEC(DOUBLE_VEC vec);

uint64_t uniform_int(uint64_t a, uint64_t b = INT_MAX);

DOUBLE_ARR mpfr_to_double_ARR(MPFR_ARR arr);

#endif
