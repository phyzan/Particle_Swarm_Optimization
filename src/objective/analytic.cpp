#include <src/objective/analytic.hpp>

#include <stdexcept>

template<typename T>
typename AnalyticObjective<T>::Function
AnalyticObjective<T>::parse(const std::string& name)
{
    if (name == "sphere")     { return Function::Sphere; }
    if (name == "rastrigin")  { return Function::Rastrigin; }
    if (name == "rosenbrock") { return Function::Rosenbrock; }

    // validate() rejects unknown names before this runs.
    throw std::runtime_error("unknown analytic function \"" + name + "\"");
}

template class AnalyticObjective<Real>;
template class AnalyticObjective<LongReal>;
