#ifndef PSO_DRIVER_SEARCH_HPP
#define PSO_DRIVER_SEARCH_HPP

#include <src/driver/dispatch.hpp>
#include <ostream>

namespace driver{

/// @brief Why a phase stopped, as text.
const char* describe(int status);

/// @brief One optimisation: fast phase, precision handoff, exact phase.
Minima run_swarm(const input::Parameters& params, std::ostream& out, Minima inherited = {});

/// @brief Runs a single swarm to convergence.
Minima run_single(const input::Parameters& params, std::ostream& out);

/// @brief Restarts repeatedly, each run repelled from the minima already found.
Minima run_deflection(const input::Parameters& params, std::ostream& out);

} // namespace driver

#endif // PSO_DRIVER_SEARCH_HPP
