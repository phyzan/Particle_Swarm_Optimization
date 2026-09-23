#ifndef PSO_DRIVER_DISPATCH_HPP
#define PSO_DRIVER_DISPATCH_HPP

#include <src/input/parameters.hpp>

#include <array>
#include <ostream>
#include <vector>

namespace driver{

// Every minimum located, in LongReal precision. One minimum occupies
// poincare.n_crossings consecutive entries.
using Minima = std::vector<std::array<LongReal, PSO_DIM>>;

/// @brief Runs the strategy the parameters select.
Minima run(const input::Parameters& params, std::ostream& out);

} // namespace driver

#endif // PSO_DRIVER_DISPATCH_HPP
