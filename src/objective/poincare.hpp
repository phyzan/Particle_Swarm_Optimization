#ifndef PSO_OBJECTIVE_POINCARE_HPP
#define PSO_OBJECTIVE_POINCARE_HPP

#include <src/integrate/integrator.hpp>
#include <src/objective/objective.hpp>

#include <array>
#include <memory>
#include <vector>

/// @brief Distance from a periodic orbit of a caldera Hamiltonian.
///
/// A particle (x, px) is lifted onto the energy shell, integrated until it
/// returns to the surface of section, and scored by how far it missed closing.
/// A zero is a periodic orbit.
///
/// The integrator is built with a null initial state, which odecraft marks as
/// not-yet-running, so a forgotten lift fails loudly rather than silently
/// integrating a meaningless state.
template<typename T>
class PoincareObjective final : public Objective<T>{

    static constexpr size_t N = 4;
    static_assert(PSO_DIM == 2, "PSO_DIM must be 2 for this particle implementation.");

public:

    // Outcome of one integration, so callers can tell "no orbit" from "error".
    enum class Status : uint8_t{
        Crossed,        // reached n_crossings
        OffSurface,     // the lifted state is not on the energy shell
        EnergyDrift,    // integration lost the energy shell
        TimedOut,       // t_max reached before n_crossings
        WrongBranch     // final py <= 0
    };

    /// @brief Builds the solver from the parameters, with no initial state yet.
    explicit PoincareObjective(const input::Parameters& pin) : 
        solver_(make_integrator<T>(pin)),
        potential{.c1=T{pin.objective.poincare.c1},
                  .c2=T{pin.objective.poincare.c2},
                  .c3=T{pin.objective.poincare.c3}},
        y_surface(pin.objective.poincare.section),
        shell_energy(pin.objective.poincare.energy),
        dt_init(pin.objective.poincare.dt_init),
        drift_tol(pin.objective.poincare.energy_drift_tol),
        t_max(pin.objective.poincare.t_max),
        n_crossings(pin.objective.poincare.n_crossings),
        metric(pin.objective.poincare.metric){}

    /// @brief Value at x, or +inf when x is infeasible.
    T evaluate(const T* x) override{
        return lift(x[0], x[1]) ? residual() : ode::inf<T>();
    }

    /// @brief Every crossing of the orbit through x.
    std::vector<std::array<T, PSO_DIM>> refine(const T* x) override{

        std::vector<std::array<T, PSO_DIM>> out;

        if (!lift(x[0], x[1]) || advance_to_crossings(&out) != Status::Crossed){
            out.clear();
        }

        return out;
    }

    /// @brief One orbit occupies n_crossings rows.
    size_t stride() const override{
        return n_crossings;
    }

private:

    /// @brief Lifts (x, px) onto the energy shell and re-arms the solver.
    bool lift(const T& x, const T& px){

        std::array<T, N> q0 = {x, y_surface, px, T{0}};

        q0[3] = 2*(shell_energy - potential.V(x, y_surface)) - px*px;

        if (q0[3] <= 0){
            return false;
        }

        q0[3] = sqrt(q0[3]);
        start = q0;

        // Reuse the step size the previous integration settled on: this
        // particle moves smoothly between iterations, so it is a good guess.
        const T guess = solver_->has_valid_ics() ? solver_->step_size() : T{dt_init};

        return solver_->set_ics(T{0}, q0.data(), guess, 1);
    }

    /// @brief Integrates to the n-th crossing and returns how far the orbit missed closing.
    T residual(){

        Status status = advance_to_crossings();

        if (status != Status::Crossed){
            return ode::inf<T>();
        }

        const T x  = last[0];
        const T px = last[2];
        const T py = last[3];

        if (py <= 0){
            return ode::inf<T>();
        }

        const T dx  = start[0] - x;
        const T dpx = start[2] - px;
        const T dpy = start[3] - py;

        if (metric == input::PoincareMetric::Euclidean){
            return sqrt(dx*dx + dpx*dpx);
        }

        return sqrt(x*x*dx*dx + px*px*dpx*dpx + py*py*dpy*dpy)
             / (x*x + px*px + py*py);
    }

    /// @brief Total energy of the solver's current state.
    T energy() const{
        const auto q = solver_->state();
        return (q[2]*q[2] + q[3]*q[3])/2 + potential.V(q[0], q[1]);
    }

    /// @brief Steps until enough crossings are seen, or the integration fails.
    Status advance_to_crossings(std::vector<std::array<T, PSO_DIM>>* record = nullptr){

        if (!solver_->is_running()){
            return Status::OffSurface;
        }

        if (record){
            record->reserve(n_crossings);
        }

        size_t seen = 0;

        while (seen < n_crossings && solver_->advance()){

            if (solver_->at_event()){

                const auto q = solver_->state();

                for (size_t i = 0; i < N; i++){
                    last[i] = q[i];
                }

                if (record){
                    record->push_back({q[0], q[2]});
                }

                seen++;
                continue;
            }

            if (abs(energy() - shell_energy) > drift_tol*abs(shell_energy)){
                return Status::EnergyDrift;
            }

            if (solver_->time() > t_max){
                return Status::TimedOut;
            }
        }

        return (seen == n_crossings) ? Status::Crossed : Status::TimedOut;
    }

    std::unique_ptr<Integrator<T>> solver_;

    // The virtual solver interface does not expose its ODE system, and the
    // energy lift in reset_to needs V(x, y). Cheaper to keep a copy than to
    // widen the interface.
    CalderaODE<T> potential;
    std::array<T, N> start{};   // lifted initial condition
    std::array<T, N> last{};    // state at the most recent crossing

    T y_surface;
    T shell_energy;
    T dt_init;
    T drift_tol;
    T t_max;

    size_t n_crossings;
    input::PoincareMetric metric;
};

// Instantiated once in poincare.cpp.
extern template class PoincareObjective<Real>;
extern template class PoincareObjective<LongReal>;

#endif // PSO_OBJECTIVE_POINCARE_HPP
