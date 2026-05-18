#include <random>

#include "functions.hpp"

// Fixed seed for reproducibility
static constexpr unsigned int FIXED_SEED = 42;
static thread_local std::mt19937 global_generator(FIXED_SEED);

// Cache the standard U(0,1) distribution - avoid recreating on every call
static thread_local std::uniform_real_distribution<double> cached_uniform_01(0.0, 1.0);

double uniform_real(double a, double b)
{
    // Return a real number using the Uniform distribution.
    // Use cached U(0,1) and scale to [a,b] to avoid distribution construction overhead.

    return a + (b - a) * cached_uniform_01(global_generator);
}

// Cache a standard normal distribution for common case
static thread_local std::normal_distribution<double> cached_standard_normal(0.0, 1.0);

double gauss_distribution(double mean, double stddev)
{
    // Return a real number using the Gaussian (Normal) distribution.
    // Use cached standard normal and transform to avoid distribution construction overhead.

    return mean + stddev * cached_standard_normal(global_generator);
}

MPFR_ARR double_to_mpfr_ARR(MPFR_ARR& mp_arr, const DOUBLE_ARR& arr)
{
    // Turn a double eigen array to an MPFR eigen array.
    mp_arr.resize(arr.rows(), arr.cols());

    for (int i = 0; i < arr.rows(); i++)
    {
        for (int j = 0; j < arr.cols(); j++)
        {
            mp_arr(i, j) = MP_REAL(arr(i, j));
        }
    }

    return mp_arr;
}

MPFR_VEC double_to_mpfr_VEC(MPFR_VEC& mp_vec, const DOUBLE_VEC& vec)
{
    // Turn a double eigen vector to an MPFR eigen vector.
    mp_vec.resize(vec.rows(), vec.cols());


    for (int i = 0; i < vec.rows(); i++)
    {
        mp_vec(i) = MP_REAL(vec(i));
    }

    return mp_vec;
}

uint64_t uniform_int(uint64_t a, uint64_t b)
{
    // Return an integer using Uniform distribution.

    std::uniform_int_distribution<uint64_t> U(a, b);

    return U(global_generator);
}

DOUBLE_ARR mpfr_to_double_ARR(DOUBLE_ARR& db_arr, const MPFR_ARR& arr)
{
    // Turn an MPFR eigen array to a double eigen array.

    db_arr.resize(arr.rows(), arr.cols());

    for (int i = 0; i < arr.rows(); i++)
    {
        for (int j = 0; j < arr.cols(); j++)
        {
            db_arr(i, j) = double(arr(i, j));
        }
    }

    return db_arr;
}
