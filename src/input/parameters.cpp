#include <src/input/parameters.hpp>

/// @brief Sets the working precision of LongReal, in bits.
void set_long_real_precision(int bits)
{
#if PSO_USE_MPREAL_SWAP
    lazex::set_default_mpreal_prec(bits);
#else
    (void)bits;   // long double has a fixed width
#endif // PSO_USE_MPREAL_SWAP
}
