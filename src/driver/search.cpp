#include <src/driver/search.hpp>
#include <src/swarm/swarm.hpp>

#include <cmath>


namespace driver{
namespace {

/// @brief Builds the objective transform that installs a pole at every known minimum.
template<typename T>
typename Swarm<T>::Transform make_deflection(const input::Deflection& cfg, const Swarm<T>& swarm)
{
    // F(x) = (f(x) + shift) * prod_i 1 / tanh(lambda * ||x - x*_i||)
    // Reshapes f, where repulsion moves x. The two compose.
    const T lambda = T{cfg.lambda};
    const T shift  = T{cfg.shift};

    return [lambda, shift, &swarm](const T& raw, size_t j) -> T {

        using std::isinf;
        using std::sqrt;
        using std::tanh;

        const auto& state = swarm.snapshot();

        if (state.found.empty() || isinf(raw)){
            return raw;
        }

        T value = raw + shift;

        for (const auto& minimum : state.found){

            T squared = T{0};

            for (size_t i = 0; i < PSO_DIM; i++){
                const T d = state.X(j, i) - minimum[i];
                squared = squared + d*d;
            }

            const T distance = sqrt(squared);

            if (distance <= 0){
                return ode::inf<T>();
            }

            value = value / tanh(lambda * distance);
        }

        return value;
    };
}

/// @brief Installs the deflection transform when the strategy calls for it.
template<typename T>
void apply_deflection(Swarm<T>& swarm, const input::Parameters& params)
{
    if (params.strategy == input::SearchStrategy::Deflection){
        swarm.set_transform(make_deflection<T>(params.deflection, swarm));
    }
}

/// @brief Narrows inherited minima to the fast phase's precision.
std::vector<std::array<Real, PSO_DIM>> narrow(const Minima& inherited)
{
    std::vector<std::array<Real, PSO_DIM>> out;
    out.reserve(inherited.size());

    for (const auto& point : inherited){
        std::array<Real, PSO_DIM> converted{};
        for (size_t i = 0; i < PSO_DIM; i++){
            converted[i] = Real(point[i]);
        }
        out.push_back(converted);
    }

    return out;
}

} // namespace

const char* describe(int status)
{
    switch (Status(status)){
    case Status::Converged:  return "converged";
    case Status::Exhausted:  return "iteration budget exhausted";
    case Status::Stalled:    return "stalled";
    case Status::Swap:       return "swapped";
    case Status::Infeasible: return "all initial particles off the energy shell";
    }
    return "unknown";
}

Minima run_swarm(const input::Parameters& params, std::ostream& out, Minima inherited)
{
    Swarm<Real> fast(params, &out);

    apply_deflection(fast, params);
    fast.initialise(narrow(inherited));

    Status status = fast.run(true);

    if (status != Status::Swap){

        if (status == Status::Converged){
            fast.record_best();
        }

        out << "  [double] " << describe(int(status))
            << " after " << fast.snapshot().iter << " iterations"
            << ", residual " << fast.residual() << '\n';

        return fast.snapshot().template cast<LongReal>().found;
    }

    out << "  [double] residual " << fast.residual()
        << " below swap_tol after " << fast.snapshot().iter
        << " iterations; continuing at " << params.convergence.mpfr_prec << " bits\n";

    Swarm<LongReal> exact(params, &out);

    apply_deflection(exact, params);

    // The handoff KEEPS the personal bests it carries across: only F is
    // recomputed at the new precision and update_bests reconciles.
    // Reassigning P = X there would discard the whole handover.
    exact.resume(fast.snapshot().template cast<LongReal>());

    status = exact.run(false);

    if (status == Status::Converged){
        exact.record_best();
    }

    out << "  [mpreal] " << describe(int(status))
        << " after " << exact.snapshot().iter << " iterations"
        << ", residual " << exact.residual() << '\n';

    return exact.snapshot().found;
}

Minima run_single(const input::Parameters& params, std::ostream& out)
{
    return run_swarm(params, out);
}

Minima run_deflection(const input::Parameters& params, std::ostream& out)
{
    Minima found;

    for (size_t run = 0; run < params.deflection.runs; run++){

        const size_t before = found.size();

        out << "-- deflection run " << (run + 1) << " of " << params.deflection.runs
            << " (" << before/objective_stride(params) << " minima so far)\n";

        found = run_swarm(params, out, std::move(found));

        if (params.deflection.stop_when_barren && found.size() == before){
            out << "-- no new minimum; stopping early\n";
            break;
        }
    }

    return found;
}

} // namespace driver
