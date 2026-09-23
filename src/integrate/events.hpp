#ifndef PSO_INTEGRATE_EVENTS_HPP
#define PSO_INTEGRATE_EVENTS_HPP

#include <src/integrate/system.hpp>

#include <string>
#include <utility>

/// @brief Event function whose zeros are crossings of the surface of section.
///
/// A concrete type rather than a std::function, for the same reason as
/// CalderaODE: the stepper checks it on every step.
template<typename T>
struct PoincareCrossing{

    T y_surface;

    /// @brief Signed distance of the state from the surface of section.
    inline T operator()(const T& /*t*/, const T* q) const{
        return q[1] - y_surface;
    }
};


/// @brief The event type the section crossing is monitored through.
template<typename T>
using SectionEvent = ode::PreciseEvent<T, PoincareCrossing<T>, ode::rhs_t<T>, ode::EventPolicy::Virtual>;

/// @brief Builds the section-crossing event the parameters describe.
template<typename T>
inline ode::EventList<T> make_section_events(const input::Parameters& pin)
{
    return ode::make_event_list<T>(
        ode::make_precise_event<T>(
            "poincare_surface",
            PoincareCrossing<T>{.y_surface = T{pin.objective.poincare.section}},
            T{pin.objective.poincare.section_ftol},
            1));
}


// One instantiation per scalar type, in events.cpp. The event machinery is
// shared by every stepper, so compiling it once keeps it off each stepper's
// critical path.
extern template struct PoincareCrossing<Real>;
extern template struct PoincareCrossing<LongReal>;

extern template ode::EventList<Real>     make_section_events<Real>(const input::Parameters&);
extern template ode::EventList<LongReal> make_section_events<LongReal>(const input::Parameters&);

#endif // PSO_INTEGRATE_EVENTS_HPP
