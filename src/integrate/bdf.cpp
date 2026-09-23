// One stepper, one translation unit. Splitting these is what lets the build
// use more than one core: each instantiates the BDF stepper for Real and for
// LongReal, with CalderaODE::Rhs inlined into it, and nothing else.

#include <src/integrate/integrator.hpp>

template class TypedIntegrator<Real, ode::Stepper::BDF>;
template class TypedIntegrator<LongReal, ode::Stepper::BDF>;
