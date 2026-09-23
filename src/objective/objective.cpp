#include <src/objective/objective.hpp>

#include <src/objective/analytic.hpp>
#include <src/objective/poincare.hpp>

/// @brief Builds the objective the parameters select.
///
/// The only place that maps an ObjectiveKind onto a concrete class. Adding an
/// objective means a new class, a new enum value, its own parameter block, and
/// one line here -- nothing in the swarm or the drivers changes.
template<typename T>
std::unique_ptr<Objective<T>> make_objective(const input::Parameters& pin)
{
    switch (pin.objective.kind){

        case input::ObjectiveKind::Analytic:
            return std::make_unique<AnalyticObjective<T>>(pin);

        case input::ObjectiveKind::Poincare:
            break;
    }

    return std::make_unique<PoincareObjective<T>>(pin);
}

/// @brief stride() without building an objective, for reporting.
size_t objective_stride(const input::Parameters& pin)
{
    switch (pin.objective.kind){
        case input::ObjectiveKind::Analytic: return 1;
        case input::ObjectiveKind::Poincare: break;
    }
    return pin.objective.poincare.n_crossings;
}

template std::unique_ptr<Objective<Real>>     make_objective<Real>(const input::Parameters&);
template std::unique_ptr<Objective<LongReal>> make_objective<LongReal>(const input::Parameters&);
