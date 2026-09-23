#include <src/input/validate.hpp>

#include <cmath>
#include <thread>

namespace input{

/// @brief Checks every cross-field rule and returns all violations.
std::vector<std::string> validate(const Parameters& p)
{
    std::vector<std::string> errors;

    const auto require = [&errors](bool ok, std::string message){
        if (!ok){
            errors.push_back(std::move(message));
        }
    };

    // ---- domain ----
    for (size_t i = 0; i < PSO_DIM; i++){
        require(p.domain.lower[i] < p.domain.upper[i],
                "domain.lower[" + std::to_string(i) + "] must be strictly less than domain.upper["
                + std::to_string(i) + "] (a zero-width box is usually a missing [domain] block)");
    }

    // ---- convergence ----
    const auto& cv = p.convergence;
    require(cv.objfun_tol > 0, "convergence.objfun_tol must be positive");
    require(cv.max_iter >= 1, "convergence.max_iter must be at least 1");
    require(cv.swap_tol > cv.objfun_tol,
            "convergence.swap_tol must exceed objfun_tol, else the run converges before it can swap");
    require(cv.swap_tol > cv.objfun_target,
            "convergence.swap_tol must exceed objfun_target, else the swap never triggers");
    require(cv.mpfr_prec >= 64, "convergence.mpfr_prec must be at least 64");
    require(std::pow(2.0, -double(cv.mpfr_prec)) < cv.objfun_tol/100,
            "convergence.mpfr_prec is too low to resolve objfun_tol");

    // ---- dynamics ----
    const auto& dy = p.dynamics;
    require(dy.pop_size > 0, "dynamics.pop_size must be positive");
    require(dy.c1 >= 0 && dy.c2 >= 0, "dynamics.c1 and c2 must be non-negative");
    require(dy.w_min >= 0, "dynamics.w_min must be non-negative");
    require(dy.w_max > dy.w_min, "dynamics.w_max must exceed w_min");
    require(dy.w_decay_frac > 0 && dy.w_decay_frac <= 1,
            "dynamics.w_decay_frac must lie in (0, 1]");
    require(dy.vclamp_frac > 0, "dynamics.vclamp_frac must be positive");

    if (dy.constriction){
        require(dy.c1 + dy.c2 > 4,
                "dynamics.constriction needs c1 + c2 > 4 (Clerc's factor is undefined otherwise)");
        require(dy.constriction_k != 0, "dynamics.constriction_k must not be zero");
    }

    if (dy.w_schedule == InertiaSchedule::Exponential){
        require(dy.w_min > 0, "dynamics.w_schedule = exponential needs w_min > 0");
    }

    // ---- topology ----
    const auto& tp = p.topology;
    require(tp.n_neighbours <= dy.pop_size,
            "topology.n_neighbours cannot exceed dynamics.pop_size");
    require(tp.radius_frac > 0 && tp.radius_frac <= 1,
            "topology.radius_frac must lie in (0, 1]");
    require(tp.dedup_tol >= 0, "topology.dedup_tol must be non-negative");

    if (tp.kind == TopologyKind::Lsh){
        require(tp.lsh.n_projections >= 4, "lsh.n_projections must be at least 4");
        require(tp.lsh.window >= 2, "lsh.window must be at least 2");
        require(tp.lsh.n_tables >= 5, "lsh.n_tables must be at least 5");
    }

    // ---- strategy ----
    if (p.strategy == SearchStrategy::Deflection){
        require(p.deflection.runs >= 1, "deflection.runs must be at least 1");
        require(p.deflection.lambda > 0, "deflection.lambda must be positive");
    }

    if (p.repulsion.enabled){
        require(p.repulsion.radius > 0, "repulsion.radius must be positive");
        require(p.repulsion.rho != 0, "repulsion.rho must not be zero");
    }

    if (p.strategy == SearchStrategy::Grid){
        require(p.grid.start_subspace < (size_t{1} << (PSO_DIM * p.grid.resolution)),
                "grid.start_subspace is beyond the last subspace");
        require(p.grid.start_subspace == 0 || p.execution.resume,
                "grid.start_subspace > 0 needs execution.resume = true, "
                "otherwise the output directory is cleared and the earlier results are lost");
    }

    // ---- objective ----
    // Only the selected objective's block is checked: the others are not read,
    // so their defaults must not be able to fail a run that never uses them.
    if (p.objective.kind == ObjectiveKind::Poincare){
        const auto& pc = p.objective.poincare;
        require(pc.n_crossings >= 1, "poincare.n_crossings must be at least 1");
        require(pc.dt_init > 0, "poincare.dt_init must be positive");
        require(pc.rtol > 0, "poincare.rtol must be positive");
        require(pc.atol >= 0, "poincare.atol must be non-negative");
        require(pc.energy_drift_tol > 0, "poincare.energy_drift_tol must be positive");
        require(pc.t_max > 0, "poincare.t_max must be positive");
        require(pc.rtol <= cv.objfun_tol/1000,
                "poincare.rtol must be at least 1000x tighter than convergence.objfun_tol, "
                "otherwise the residual is integrator noise at the scale being resolved");

        const bool known_integrator =
            pc.integrator == "Euler"  || pc.integrator == "RK4"    || pc.integrator == "RK23" ||
            pc.integrator == "RK45"   || pc.integrator == "DOP853" || pc.integrator == "BDF";

        require(known_integrator,
                "poincare.integrator is \"" + pc.integrator
                + "\"; expected Euler, RK4, RK23, RK45, DOP853 or BDF");

        if (p.grid.batch > 0 && p.strategy == SearchStrategy::Grid){
            const size_t hardware = std::max<size_t>(1, std::thread::hardware_concurrency());
            require(p.grid.batch * p.execution.threads <= hardware,
                    "grid.batch * execution.threads exceeds the hardware thread count; "
                    "the two multiply, so the machine ends up oversubscribed "
                    "(leave grid.batch at 0 to have it derived)");
        }
    }

    if (p.objective.kind == ObjectiveKind::Analytic){
        const std::string& fn = p.objective.analytic.function;
        require(fn == "sphere" || fn == "rastrigin" || fn == "rosenbrock",
                "analytic.function is \"" + fn
                + "\"; expected sphere, rastrigin or rosenbrock");
    }

    // ---- execution ----
    require(p.execution.threads >= 1, "execution.threads must be at least 1");

    return errors;
}

} // namespace input
