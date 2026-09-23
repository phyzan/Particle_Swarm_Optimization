#ifndef PSO_INTEGRATE_INTEGRATOR_HPP
#define PSO_INTEGRATE_INTEGRATOR_HPP

#include <src/integrate/events.hpp>
#include <src/integrate/system.hpp>

#include <memory>

/// @brief What a particle needs of an integrator, with the stepper erased.
///
/// ode's steppers take the method as a template parameter, so a concrete
/// class per stepper is unavoidable. This interface keeps that fact local:
/// one virtual call per particle operation -- not per ODE step -- and
/// CalderaODE::Rhs still inlines inside the stepper.
template<typename T>
class Integrator{
public:

    virtual ~Integrator() = default;

    /// @brief Re-arms the solver from a new initial condition.
    virtual bool set_ics(const T& t0, const T* q, const T& stepsize, int dir) = 0;

    /// @brief Takes one step; false when the solver can go no further.
    virtual bool advance() = 0;

    /// @brief Whether the step just taken landed on the section crossing.
    virtual bool at_event() const = 0;

    /// @brief The current state vector.
    virtual const T* state() const = 0;

    /// @brief The current time.
    virtual const T& time() const = 0;

    /// @brief Whether the solver holds usable initial conditions.
    virtual bool has_valid_ics() const = 0;

    /// @brief Whether the solver is alive and able to advance.
    virtual bool is_running() const = 0;

    /// @brief The step size the last integration settled on.
    virtual T step_size() const = 0;
};


/// @brief One concrete stepper, wrapped behind Integrator.
template<typename T, ode::Stepper S>
class TypedIntegrator final : public Integrator<T>{

    static constexpr size_t N = 4;

    using Solver = typename ode::detail::SolverTypeGetter<
        S, T, N, ode::SolverPolicy::RichStatic, CalderaODE<T>, void>::type;

public:

    /// @brief Builds the stepper with no initial state yet.
    explicit TypedIntegrator(const input::Parameters& pin)
        : solver_(make_system<T>(pin),
                  T{0},                                   // t0
                  ode::View1D<T, N>{nullptr},             // q0, set by set_ics
                  T{pin.objective.poincare.rtol},
                  T{pin.objective.poincare.atol},
                  T{0},                                   // min_step
                  T{0},                                   // max_step (0 -> inf)
                  T{pin.objective.poincare.dt_init},
                  1,                                      // direction
                  make_section_events<T>(pin)){}

    bool set_ics(const T& t0, const T* q, const T& stepsize, int dir) override{
        return solver_.set_ics(t0, q, stepsize, dir);
    }

    bool advance() override{
        return solver_.advance();
    }

    bool at_event() const override{
        return solver_.at_event();
    }

    const T* state() const override{
        return solver_.vector().data();
    }

    const T& time() const override{
        return solver_.t();
    }

    bool has_valid_ics() const override{
        return solver_.has_valid_ics();
    }

    bool is_running() const override{
        return solver_.is_running();
    }

    T step_size() const override{
        return solver_.ics().habs();
    }

private:
    Solver solver_;
};


/// @brief Builds the integrator for the stepper named in the parameters.
template<typename T>
std::unique_ptr<Integrator<T>> make_integrator(const input::Parameters& pin);


// Each of these is instantiated in its own translation unit under integrate/,
// so the twelve stepper instantiations compile in parallel rather than
// piling into whichever file calls make_integrator.
#define PSO_EXTERN_STEPPER(S)                                   \
    extern template class TypedIntegrator<Real, S>;             \
    extern template class TypedIntegrator<LongReal, S>;

PSO_EXTERN_STEPPER(ode::Stepper::Euler)
PSO_EXTERN_STEPPER(ode::Stepper::RK4)
PSO_EXTERN_STEPPER(ode::Stepper::RK23)
PSO_EXTERN_STEPPER(ode::Stepper::RK45)
PSO_EXTERN_STEPPER(ode::Stepper::DOP853)
PSO_EXTERN_STEPPER(ode::Stepper::BDF)

#undef PSO_EXTERN_STEPPER

extern template std::unique_ptr<Integrator<Real>>     make_integrator<Real>(const input::Parameters&);
extern template std::unique_ptr<Integrator<LongReal>> make_integrator<LongReal>(const input::Parameters&);

#endif // PSO_INTEGRATE_INTEGRATOR_HPP
