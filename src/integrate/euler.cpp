// One stepper, one translation unit. Splitting these is what lets the build
// use more than one core: each instantiates the Euler stepper for Real and for
// LongReal, with CalderaODE::Rhs inlined into it, and nothing else.

#include <src/integrate/integrator.hpp>

template class TypedIntegrator<Real, ode::Stepper::Euler>;
template class TypedIntegrator<LongReal, ode::Stepper::Euler>;
