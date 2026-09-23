#include <src/swarm/swarm.hpp>

// The one place Swarm is instantiated. Member templates -- SwarmState::cast
// among them -- are not covered by an explicit class instantiation, so those
// are still instantiated where they are used.
template struct SwarmState<Real>;
template struct SwarmState<LongReal>;
template class Swarm<Real>;
template class Swarm<LongReal>;
