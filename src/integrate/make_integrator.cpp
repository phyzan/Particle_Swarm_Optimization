#include <src/integrate/integrator.hpp>

#include <stdexcept>

/// @brief Builds the integrator for the stepper named in the parameters.
///
/// choose_integrator_case is odecraft's runtime-to-compile-time switch, and
/// this is the single place in the program where that crossing happens. The
/// extern template declarations in the header keep it from instantiating any
/// stepper here: each one is already compiled in its own translation unit.
template<typename T>
std::unique_ptr<Integrator<T>> make_integrator(const input::Parameters& pin)
{
    return ode::choose_integrator_case<std::unique_ptr<Integrator<T>>>(
        ode::getIntegrator(pin.objective.poincare.integrator),
        [&pin]<ode::Stepper S>(){
            return std::unique_ptr<Integrator<T>>(new TypedIntegrator<T, S>(pin));
        });
}

template std::unique_ptr<Integrator<Real>>     make_integrator<Real>(const input::Parameters&);
template std::unique_ptr<Integrator<LongReal>> make_integrator<LongReal>(const input::Parameters&);
