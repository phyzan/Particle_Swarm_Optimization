#include <src/integrate/events.hpp>

template struct PoincareCrossing<Real>;
template struct PoincareCrossing<LongReal>;

template ode::EventList<Real>     make_section_events<Real>(const input::Parameters&);
template ode::EventList<LongReal> make_section_events<LongReal>(const input::Parameters&);
