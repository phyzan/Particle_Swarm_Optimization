#ifndef PSO_DRIVER_GRID_HPP
#define PSO_DRIVER_GRID_HPP

#include <src/driver/dispatch.hpp>

#include <array>
#include <cstdint>
#include <ostream>


namespace driver{

/// @brief Number of subspaces the resolution produces.
size_t grid_size(const input::Grid& cfg);

/// @brief Bounds of one subspace, by mixed-radix decomposition of its index.
void subspace_bounds(size_t cell, const input::Parameters& params,
                     std::array<Real, PSO_DIM>& lower,
                     std::array<Real, PSO_DIM>& upper);

/// @brief Derives an independent RNG stream for one worker.
uint64_t seed_for(uint64_t master, uint64_t worker);

/// @brief Runs one swarm per subspace, `batch` at a time, and merges the results.
Minima run_grid(const input::Parameters& params, std::ostream& out);

} // namespace driver

#endif // PSO_DRIVER_GRID_HPP
