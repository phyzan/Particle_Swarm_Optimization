#include <random>

#include "functions.hpp"

double uniform_real(double a, double b)
{
    // Return a real number using the Uniform distribution.

    std::random_device random_seed;
    std::default_random_engine generator(random_seed());
    std::uniform_real_distribution<double> U(a, b);

    return U(generator);
}

double gauss_distribution(double mean, double stddev)
{
    // Return a real number using the Gaussian (Normal) distribution.

    std::random_device random_seed;
    std::mt19937 generator(random_seed());
    std::normal_distribution<double> Gauss(mean, stddev);

    return Gauss(generator);
}

MPFR_ARR double_to_mpfr_ARR(DOUBLE_ARR arr)
{
    // Turn a double eigen array to an MPFR eigen array.

    MPFR_ARR mp_arr(arr.rows(), arr.cols());

    for (int i = 0; i < arr.rows(); i++)
    {
        for (int j = 0; j < arr.cols(); j++)
        {
            mp_arr(i, j) = MP_REAL(arr(i, j));
        }
    }

    return mp_arr;
}

MPFR_VEC double_to_mpfr_VEC(DOUBLE_VEC vec)
{
    // Turn a double eigen vector to an MPFR eigen vector.

    MPFR_VEC mp_vec(vec.rows(), vec.cols());

    for (int i = 0; i < vec.rows(); i++)
    {
        mp_vec(i) = MP_REAL(vec(i));
    }

    return mp_vec;
}

uint64_t uniform_int(uint64_t a, uint64_t b)
{
    // Return an integer using Uniform distribution.

    std::random_device random_seed;
    std::default_random_engine generator(random_seed());
    std::uniform_int_distribution<uint64_t> U(a, b);

    return U(generator);
}

DOUBLE_ARR mpfr_to_double_ARR(MPFR_ARR arr)
{
    // Turn an MPFR eigen array to a double eigen array.

    DOUBLE_ARR db_arr(arr.rows(), arr.cols());

    for (int i = 0; i < arr.rows(); i++)
    {
        for (int j = 0; j < arr.cols(); j++)
        {
            db_arr(i, j) = double(arr(i, j));
        }
    }

    return db_arr;
}
