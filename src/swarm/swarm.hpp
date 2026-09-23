#ifndef PSO_SWARM_SWARM_HPP
#define PSO_SWARM_SWARM_HPP

#include <src/input/parameters.hpp>
#include <src/objective/objective.hpp>
#include <src/swarm/topology.hpp>

#include <array>
#include <memory>
#include <cmath>
#include <odecraft/odecraft.hpp>
#include <functional>
#include <random>
#include <vector>

#if defined(PSO_USE_GCD)
#include <dispatch/dispatch.h>
#elif defined(_OPENMP)
#include <omp.h>
#endif

using ode::Array1D;
using ode::Array2D;


// Why a phase stopped. SWAP is not a failure: it means the residual dropped
// below swap_tol and the run should continue at higher precision.
enum class Status : uint8_t{
    Converged,
    Exhausted,
    Stalled,
    Swap,
    Infeasible
};


// Everything that must survive the precision swap, and nothing else.
// Particles are ROWS: (j, i) is dimension i of particle j. ndspan defaults to
// row-major, so one particle's coordinates sit contiguously.
template<typename T>
struct SwarmState{

    Array2D<T, 0, PSO_DIM> X;    // current positions
    Array2D<T, 0, PSO_DIM> V;    // current velocities
    Array2D<T, 0, PSO_DIM> P;    // personal-best positions
    Array1D<T> F;                // f(X[j])
    Array1D<T> Fp;               // f(P[j]) -- copied from F on improvement

    // Minima accumulated so far. Each one occupies n_crossings consecutive
    // entries, so the count of minima is size()/n_crossings.
    std::vector<std::array<T, PSO_DIM>> found;

    size_t g = 0;       // argmin(Fp)
    size_t iter = 0;    // carried across the swap: max_iter is a shared budget
    size_t stall = 0;   // iterations since Fp last improved

    /// @brief Allocates every array for a population of the given size.
    explicit SwarmState(size_t pop_size = 0)
        : X(pop_size, PSO_DIM), V(pop_size, PSO_DIM), P(pop_size, PSO_DIM),
          F(pop_size), Fp(pop_size) {}

    /// @brief Number of particles in the swarm.
    size_t pop_size() const{
        return X.Nrows();
    }

    /// @brief Converts the whole state to another precision, mapping +inf to +inf.
        template<typename U>
    SwarmState<U> cast() const{

        SwarmState<U> out(pop_size());

        for (size_t j = 0; j < pop_size(); j++){
            for (size_t i = 0; i < PSO_DIM; i++){
                out.X(j, i) = convert<U>(X(j, i));
                out.V(j, i) = convert<U>(V(j, i));
                out.P(j, i) = convert<U>(P(j, i));
            }
            out.F[j]  = convert<U>(F[j]);
            out.Fp[j] = convert<U>(Fp[j]);
        }

        out.found.reserve(found.size());

        for (const auto& point : found){
            std::array<U, PSO_DIM> converted{};
            for (size_t i = 0; i < PSO_DIM; i++){
                converted[i] = convert<U>(point[i]);
            }
            out.found.push_back(converted);
        }

        out.g = g;
        out.iter = iter;
        out.stall = stall;

        return out;
    }

private:

    /// @brief Converts one scalar, preserving infinity exactly.
    template<typename U, typename V_>
    static U convert(const V_& value){
        // std two-step: the using-declaration supplies the arithmetic overloads,
        // ADL supplies mpfr::isinf for the mpreal wrappers. An unqualified call
        // alone finds neither for double, because <cmath> only guarantees std::.
        using std::isinf;
        return isinf(value) ? ode::inf<U>() : U(value);
    }
};


// One swarm at one precision. Knows nothing about strategies or about the
// other precision -- the driver owns both of those.
template<typename T>
class Swarm{
public:

    /// @brief A transform applied to every objective value, as deflection uses.
    using Transform = std::function<T(const T& raw, size_t j)>;

    /// @brief Builds the swarm and one solver per particle.
    Swarm(const input::Parameters& pin, std::ostream* stream)
        : params(pin),
          out(stream),
          rng(pin.execution.seed),
          topo_rng(pin.execution.seed ^ 0x5DEECE66DULL),
          state(pin.dynamics.pop_size),
          neighbours(pin.topology, pin.dynamics.pop_size, span_of(pin))
    {
        for (size_t i = 0; i < PSO_DIM; i++){
            span[i]  = pin.domain.upper[i] - pin.domain.lower[i];
            v_max[i] = T{pin.dynamics.vclamp_frac * span[i]};
        }

        objectives.reserve(pop_size());

        for (size_t j = 0; j < pop_size(); j++){
            objectives.push_back(make_objective<T>(pin));
        }
    }

    /// @brief Draws a fresh population, keeping any inherited minima.
    void initialise(std::vector<std::array<T, PSO_DIM>> inherited = {}){

        std::uniform_real_distribution<Real> uniform(0.0, 1.0);

        for (size_t j = 0; j < pop_size(); j++){
            for (size_t i = 0; i < PSO_DIM; i++){

                state.X(j, i) = T{params.domain.lower[i] + uniform(rng)*span[i]};

                switch (params.dynamics.velocity_init){
                case input::VelocityInit::Zero:
                    state.V(j, i) = T{0};
                    break;
                case input::VelocityInit::ScaledUniform:
                    state.V(j, i) = T{2*uniform(rng) - 1} * v_max[i];
                    break;
                case input::VelocityInit::Uniform01:
                    state.V(j, i) = T{uniform(rng)};
                    break;
                }
            }
        }

        state.found = std::move(inherited);
        state.P = state.X;

        for (size_t j = 0; j < pop_size(); j++){
            state.Fp[j] = ode::inf<T>();
        }

        state.g = 0;
        state.iter = 0;
        state.stall = 0;

        repulsion();
        evaluate();

        state.P = state.X;
        state.Fp = state.F;

        select_best();
    }

    /// @brief Resumes from a state handed over by another precision.
    void resume(SwarmState<T> handover){

        // The personal bests carried across are KEPT: only F is recomputed
        // below, and update_bests reconciles. Overwriting P/Fp here would
        // discard the entire handover.
        state = std::move(handover);

        evaluate();
        update_bests();
    }

    /// @brief Iterates until it converges, runs out of budget, stalls, or swaps.
    Status run(bool allow_swap){

        if (!feasible()){
            return Status::Infeasible;
        }

        while (true){

            if (converged())              { return Status::Converged; }
            if (exhausted())              { return Status::Exhausted; }
            if (stalled())                { return Status::Stalled; }
            if (allow_swap && swapping()) { return Status::Swap; }

            state.iter++;

            velocity_update(inertia());
            position_update();
            evaluate();
            update_bests();
        }
    }

    /// @brief Records the best particle's full orbit as a new minimum.
    void record_best(){

        for (const auto& point : objectives[state.g]->refine(&state.P(state.g, 0))){
            state.found.push_back(point);
        }
    }

    /// @brief Installs a transform applied to every objective value.
    void set_transform(Transform t){
        transform = std::move(t);
    }

    /// @brief Distance of the best particle's value from the target.
    T residual() const{
        // std two-step. Unqualified abs() on a double picks int abs(int) from
        // <cstdlib> and truncates the residual to 0, which reads as instant
        // convergence. ADL still supplies mpfr::abs for the mpreal wrappers.
        using std::abs;
        return abs(state.Fp[state.g] - T{params.convergence.objfun_target});
    }

    /// @brief Whether the residual has reached objfun_tol.
    bool converged() const{ return residual() <= T{params.convergence.objfun_tol}; }
    /// @brief Whether the iteration budget is spent.
    bool exhausted() const{ return state.iter  >= params.convergence.max_iter; }
    /// @brief Whether the best value has not improved for max_stall_iter iterations.
    bool stalled()   const{ return state.stall >= params.convergence.max_stall_iter; }
    /// @brief Whether the residual has dropped below swap_tol.
    bool swapping()  const{ return residual() < T{params.convergence.swap_tol}; }

    /// @brief Number of particles in the swarm.
    size_t pop_size() const{ return params.dynamics.pop_size; }

    /// @brief Read-only access to the swarm state.
    const SwarmState<T>& snapshot() const{ return state; }

    /// @brief Mutable access to the swarm state.
    SwarmState<T>& snapshot(){ return state; }

    /// @brief The stream this swarm reports to.
    std::ostream& log() const{ return *out; }

private:

    /// @brief Width of the domain box in each dimension.
    static std::array<Real, PSO_DIM> span_of(const input::Parameters& pin){
        std::array<Real, PSO_DIM> out{};
        for (size_t i = 0; i < PSO_DIM; i++){
            out[i] = pin.domain.upper[i] - pin.domain.lower[i];
        }
        return out;
    }

    /// @brief Whether any particle reached the energy shell.
    bool feasible() const{
        using std::isinf;
        for (size_t j = 0; j < pop_size(); j++){
            if (!isinf(state.F[j])){
                return true;
            }
        }
        return false;
    }

    /// @brief Evaluates the objective once per particle; the only place it is called.
    void evaluate(){

        // Particle idx owns objectives[idx] and writes only state.F[idx], so
        // the loop partitioning is the whole synchronisation story: no two
        // workers ever touch the same object.
        const auto evaluate_one = [this](size_t idx){
            T value = objectives[idx]->evaluate(&state.X(idx, 0));

            // CONTRACT: every slot written, on every path.
            state.F[idx] = transform ? transform(value, idx) : value;
        };

#if defined(PSO_USE_GCD)

        // Apple clang ships no OpenMP runtime, but libdispatch is always
        // there. dispatch_apply is a parallel for and blocks until every
        // iteration is done, so it needs no group, no completion handler and
        // no per-particle heap allocation.
        //
        // Concurrency is the system's to decide here: libdispatch sizes the
        // pool from the machine and the current load, so execution.threads is
        // advisory on this path rather than a cap.
        dispatch_apply(pop_size(), DISPATCH_APPLY_AUTO, ^(size_t idx){
            evaluate_one(idx);
        });

#elif defined(_OPENMP)

        const int threads = int(params.execution.threads);

        #pragma omp parallel for schedule(dynamic) num_threads(threads)
        for (long long j = 0; j < (long long)pop_size(); j++){
            evaluate_one(size_t(j));
        }

#else

        for (size_t idx = 0; idx < pop_size(); idx++){
            evaluate_one(idx);
        }

#endif
    }

    /// @brief Adopts improved positions; the only writer of P, Fp, g and stall.
    void update_bests(){

        // Being the only writer is what gives the invariant
        // Fp[j] == f(P[j]) exactly one place to be proved.
        bool improved = false;

        for (size_t j = 0; j < pop_size(); j++){
            if (state.F[j] < state.Fp[j]){
                state.Fp[j] = state.F[j];
                for (size_t i = 0; i < PSO_DIM; i++){
                    state.P(j, i) = state.X(j, i);
                }
                improved = true;
            }
        }

        select_best();

        state.stall = improved ? 0 : state.stall + 1;
    }

    /// @brief Recomputes g as the index of the lowest personal best.
    void select_best(){
        size_t best = 0;
        for (size_t j = 1; j < pop_size(); j++){
            if (state.Fp[j] < state.Fp[best]){
                best = j;
            }
        }
        state.g = best;
    }

    /// @brief Inertia weight for the current iteration, per the schedule.
    T inertia() const{

        const auto& dyn = params.dynamics;

        if (dyn.w_schedule == input::InertiaSchedule::Constant){
            return T{dyn.w_max};
        }

        const Real decay = std::max<Real>(1, std::round(dyn.w_decay_frac * Real(params.convergence.max_iter)));
        const Real s = std::min<Real>(Real(state.iter), decay) / decay;

        if (dyn.w_schedule == input::InertiaSchedule::Exponential && dyn.w_min > 0){
            return T{dyn.w_min * std::pow(dyn.w_max / dyn.w_min, 1 - s)};
        }

        return T{dyn.w_max - s*(dyn.w_max - dyn.w_min)};
    }

    /// @brief Index of the particle supplying j's social attractor.
    size_t attractor(size_t j){

        if (neighbours.mode() == input::TopologyKind::Global){
            return state.g;
        }

        neighbours.query(state.X, j, pop_size(), candidates);

        if (candidates.empty()){
            return state.g;
        }

        // The neighbourhood is {j} together with its candidates. Seeding the
        // search with state.g instead would be unsatisfiable: g is argmin(Fp)
        // by definition, so no q could ever beat it and every topology would
        // silently collapse into Global.
        size_t best = j;
        T best_value = state.Fp[j];

        for (size_t q : candidates){
            if (state.Fp[q] < best_value){
                best_value = state.Fp[q];
                best = q;
            }
        }

        return best;
    }

    /// @brief Applies the PSO recurrence, constriction and the speed cap.
    void velocity_update(const T& w){

        const auto& dyn = params.dynamics;

        T chi = T{1};

        if (dyn.constriction){
            const Real phi = dyn.c1 + dyn.c2;
            chi = T{std::abs(2*dyn.constriction_k / (2 - phi - std::sqrt(phi*phi - 4*phi)))};
        }

        neighbours.rebuild(state.X, pop_size(), topo_rng);

        std::uniform_real_distribution<Real> uniform(0.0, 1.0);

        for (size_t j = 0; j < pop_size(); j++){

            const size_t a = attractor(j);

            for (size_t i = 0; i < PSO_DIM; i++){

                // Fresh per component: one scalar per particle would confine
                // the pull to the line joining it to its attractors.
                const T r1 = T{uniform(rng)};
                const T r2 = T{uniform(rng)};

                state.V(j, i) = w * state.V(j, i)
                              + T{dyn.c1} * r1 * (state.P(j, i) - state.X(j, i))
                              + T{dyn.c2} * r2 * (state.P(a, i) - state.X(j, i));

                if (dyn.constriction){
                    state.V(j, i) = chi * state.V(j, i);
                }

                if (state.V(j, i) < -v_max[i]){ state.V(j, i) = -v_max[i]; }
                if (state.V(j, i) >  v_max[i]){ state.V(j, i) =  v_max[i]; }
            }
        }
    }

    /// @brief Advances positions, applies the bound policy, then repulsion.
    void position_update(){

        std::uniform_real_distribution<Real> uniform(0.0, 1.0);

        for (size_t j = 0; j < pop_size(); j++){
            for (size_t i = 0; i < PSO_DIM; i++){

                state.X(j, i) = state.X(j, i) + state.V(j, i);

                const T low  = T{params.domain.lower[i]};
                const T high = T{params.domain.upper[i]};

                if (state.X(j, i) >= low && state.X(j, i) <= high){
                    continue;
                }

                switch (params.domain.bound_policy){

                case input::BoundPolicy::None:
                    break;

                case input::BoundPolicy::Clamp:
                    state.X(j, i) = (state.X(j, i) < low) ? low : high;
                    state.V(j, i) = T{0};
                    break;

                case input::BoundPolicy::Reflect:
                    state.X(j, i) = (state.X(j, i) < low) ? 2*low - state.X(j, i)
                                                          : 2*high - state.X(j, i);
                    if (state.X(j, i) < low || state.X(j, i) > high){
                        state.X(j, i) = (state.X(j, i) < low) ? low : high;
                    }
                    state.V(j, i) = -state.V(j, i);
                    break;

                case input::BoundPolicy::ReInit:
                    state.X(j, i) = T{params.domain.lower[i] + uniform(rng)*span[i]};
                    state.V(j, i) = T{0};
                    break;
                }
            }
        }

        repulsion();
    }

    /// @brief Pushes particles out of the repulsion radius of known minima.
    void repulsion(){

        if (!params.repulsion.enabled || state.found.empty()){
            return;
        }

        const T radius = T{params.repulsion.radius};
        const T rho    = T{params.repulsion.rho};

        for (const auto& minimum : state.found){
            for (size_t j = 0; j < pop_size(); j++){

                T squared = T{0};

                for (size_t i = 0; i < PSO_DIM; i++){
                    const T d = state.X(j, i) - minimum[i];
                    squared = squared + d*d;
                }

                using std::sqrt;
                const T distance = sqrt(squared);

                if (distance <= 0 || distance > radius){
                    continue;
                }

                for (size_t i = 0; i < PSO_DIM; i++){
                    state.X(j, i) = state.X(j, i) + rho*(state.X(j, i) - minimum[i])/distance;
                    state.P(j, i) = state.X(j, i);
                }

                // The particle was moved, so its recorded best no longer
                // describes anywhere it has been.
                state.Fp[j] = ode::inf<T>();
            }
        }
    }

    input::Parameters params;
    std::ostream* out;
    std::mt19937_64 rng;       // swarm dynamics
    std::mt19937_64 topo_rng;  // neighbourhood construction, kept separate

    SwarmState<T> state;

    std::array<Real, PSO_DIM> span{};
    std::array<T, PSO_DIM> v_max{};

    topology::Neighbourhood neighbours;
    std::vector<size_t> candidates;

    // One objective instance per swarm particle. An implementation may hold
    // heavy mutable state -- the Poincare one owns an ODE solver -- and the
    // parallel-for partitions j across threads, so no two threads ever touch
    // the same instance.
    std::vector<std::unique_ptr<Objective<T>>> objectives;

    Transform transform;
};

// Instantiated once in swarm.cpp. Without these, every translation unit that
// touches a Swarm instantiates the whole class again -- which is most of what
// the build spends its time on.
extern template struct SwarmState<Real>;
extern template struct SwarmState<LongReal>;
extern template class Swarm<Real>;
extern template class Swarm<LongReal>;

#endif // PSO_SWARM_SWARM_HPP
