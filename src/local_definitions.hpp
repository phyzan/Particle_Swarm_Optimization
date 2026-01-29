#ifndef __MP_DEF__
#define __MP_DEF__

#include <cmath>
#include <limits>
#include <unsupported/Eigen/MPRealSupport>

#ifndef MP_REAL
#define MP_REAL mpfr::mpreal
#endif

#ifndef MPFR_VEC
#define MPFR_VEC Eigen::Vector<MP_REAL, Eigen::Dynamic>
#endif

#ifndef MPFR_ARR
#define MPFR_ARR Eigen::Array<MP_REAL, Eigen::Dynamic, Eigen::Dynamic>
#endif

#ifndef MPFR_EMPTY
#define MPFR_EMPTY Eigen::Array<MP_REAL, Eigen::Dynamic, 1>
#endif

#ifndef DOUBLE_VEC
#define DOUBLE_VEC Eigen::Vector<double, Eigen::Dynamic>
#endif

#ifndef DOUBLE_ARR
#define DOUBLE_ARR Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>
#endif

#ifndef DOUBLE_EMPTY
#define DOUBLE_EMPTY Eigen::Array<double, Eigen::Dynamic, 1>
#endif

namespace lmath
{ // Local custom math functions.

inline double abs(double x)
{
    return std::abs(x);
}

inline double sqrt(double x)
{
    return std::sqrt(x);
}

inline double floor(double x)
{
    return std::floor(x);
}

inline double pow(double base, double exp)
{
    return std::pow(base, exp);
}

inline mpfr::mpreal abs(const mpfr::mpreal &x)
{
    return mpfr::abs(x);
}

inline mpfr::mpreal sqrt(const mpfr::mpreal &x)
{
    return mpfr::sqrt(x);
}

inline mpfr::mpreal floor(const mpfr::mpreal &x)
{
    return mpfr::floor(x);
}

inline mpfr::mpreal pow(const mpfr::mpreal &base, const mpfr::mpreal &exp)
{
    return mpfr::pow(base, exp);
}

template <typename Type> inline Type get_infinity();

template <> inline double get_infinity<double>()
{
    return std::numeric_limits<double>::infinity();
}

template <> inline mpfr::mpreal get_infinity<mpfr::mpreal>()
{
    return mpfr::mpreal("inf");
}

template <typename Type> bool isinf(const Type &x);

template <> inline bool isinf<double>(const double &x)
{
    return std::isinf(x);
}

template <> inline bool isinf<mpfr::mpreal>(const mpfr::mpreal &x)
{
    return mpfr::isinf(x) != 0;
}

} // namespace lmath

#endif