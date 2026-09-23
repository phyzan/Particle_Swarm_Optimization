#ifndef PSO_OBJECTIVE_OBJECTIVE_HPP
#define PSO_OBJECTIVE_OBJECTIVE_HPP

#include <src/input/parameters.hpp>

#include <odecraft/odecraft.hpp>

#include <array>
#include <memory>
#include <vector>

/// @brief Everything the swarm needs of an objective function.
///
/// The swarm touches an objective in exactly two ways: it asks for a value at
/// a position, and -- once, when a minimum is found -- for the full record of
/// the solution there. Nothing above this interface knows what is being
/// minimised, which is why the drivers, the precision swap, deflection,
/// repulsion and the grid are all objective-agnostic.
///
/// One instance per swarm particle: an implementation may carry heavy mutable
/// state (the Poincare one owns an ODE solver), and the evaluation loop runs
/// particles in parallel.
template<typename T>
class Objective{
public:

    virtual ~Objective() = default;

    /// @brief Value at x, or +inf when x is infeasible.
    ///
    /// Called pop_size times per iteration -- this is the hot path, and on
    /// most problems it dominates everything else the program does.
    virtual T evaluate(const T* x) = 0;

    /// @brief The full solution at x, recorded once a minimum is accepted.
    ///
    /// Returns stride() rows, or nothing when x does not in fact solve the
    /// problem. Called once per minimum, never in the evaluation loop, so it
    /// should be written for clarity rather than speed.
    virtual std::vector<std::array<T, PSO_DIM>> refine(const T* x) = 0;

    /// @brief How many rows one recorded minimum occupies.
    virtual size_t stride() const = 0;
};


/// @brief Builds the objective the parameters select.
template<typename T>
std::unique_ptr<Objective<T>> make_objective(const input::Parameters& pin);

/// @brief stride() without building an objective, for reporting.
size_t objective_stride(const input::Parameters& pin);

extern template std::unique_ptr<Objective<Real>>     make_objective<Real>(const input::Parameters&);
extern template std::unique_ptr<Objective<LongReal>> make_objective<LongReal>(const input::Parameters&);

#endif // PSO_OBJECTIVE_OBJECTIVE_HPP
