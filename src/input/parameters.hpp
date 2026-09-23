#ifndef PSO_INPUT_PARAMETERS_HPP
#define PSO_INPUT_PARAMETERS_HPP

#include <array>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <string>

#ifndef PSO_DIM
#define PSO_DIM 2
#endif // PSO_DIM

#ifndef PSO_USE_MPREAL_SWAP
#define PSO_USE_MPREAL_SWAP true
#endif // PSO_USE_MPREAL_SWAP

using Real = double;

#if PSO_USE_MPREAL_SWAP
#include <lazex/apps/lazex_mpreal.hpp>
using LongReal = lazex::LazyType<mpfr::mpreal>;
#else
using LongReal = long double;
#endif // PSO_USE_MPREAL_SWAP

/// @brief Sets the working precision of LongReal, in bits.
void set_long_real_precision(int bits);

namespace input{

// ---------------------------------------------------------------------------
// Search domain
// ---------------------------------------------------------------------------

// What happens to a particle that has left the search box after a position
// update. The box always seeds the initial population and sets the velocity
// clamp; this decides whether it also confines the swarm during the run.
enum class BoundPolicy : uint8_t{
    None,
    Clamp,
    Reflect,
    ReInit
};

struct Domain{
    std::array<Real, PSO_DIM> lower; // Constraint: lower < upper
    std::array<Real, PSO_DIM> upper;
    BoundPolicy bound_policy = BoundPolicy::None;
    // None reproduces the reference implementation: the box seeds the initial
    // population and sets the velocity clamp, but particles may leave it. The
    // objective returning +inf outside the feasible set acts as a soft wall.
};

// ---------------------------------------------------------------------------
// Termination
// ---------------------------------------------------------------------------

struct SwarmConvergence{
    Real objfun_target = 0.0;
    // The objective VALUE to reach, not its location. For a return-map residual
    // a periodic orbit has residual exactly 0, so this is a genuine known input.
    Real objfun_tol = 1e-8; // Constraint: tol > 0
    // Success test: |f_best - objfun_target| <= objfun_tol.

    Real swap_tol = 1e-5; // Constraint: objfun_target < swap_tol < start residual
    // Precision handoff point. The run begins in Real and switches to LongReal
    // once |f_best - objfun_target| < swap_tol, carrying positions, velocities,
    // personal bests and the iteration counter across untouched. Real is fast
    // but its own round-off swamps the residual below ~1e-6; LongReal reaches
    // objfun_tol but costs far more per evaluation (~200x for mpreal). Degenerate
    // settings are useful: <= objfun_target never swaps (pure Real run), while a
    // value above the starting residual runs entirely in LongReal.
    // Constraint: swap_tol > objfun_tol, else the run converges before swapping.

    int mpfr_prec = 100; // Constraint: >= 64, and 2^-mpfr_prec << objfun_tol
    // Working precision in bits for LongReal, applied once via
    // lazex::set_default_mpreal_prec(mpfr_prec) before any LongReal is
    // constructed. Ignored when PSO_USE_MPREAL_SWAP is off, since long double
    // has a fixed width. The default leaves ~8e-31 of headroom under the 1e-8
    // objfun_tol; raising it costs runtime roughly linearly.

    size_t max_iter = 15000; // Constraint: max_iter >= 1
    // Shared budget across both precision phases, not per phase.
    size_t max_stall_iter = std::numeric_limits<size_t>::max(); // Constraint: >= 0
};

// ---------------------------------------------------------------------------
// Swarm dynamics
// ---------------------------------------------------------------------------

// How the inertia weight w is driven from w_max down to w_min over a run.
// w sets the exploration/exploitation balance: high w lets particles coast on
// their own momentum, low w lets the c1/c2 attractors dominate.
enum class InertiaSchedule : uint8_t{
    Linear,
    Exponential,
    Constant
};

// How the velocity field is drawn at iteration 0. Only shapes the first few
// iterations, but a scale-blind choice biases early motion when the domain's
// dimensions have very different spans.
enum class VelocityInit : uint8_t{
    Zero,
    ScaledUniform, // U(-v_max, +v_max), per dimension
    Uniform01      // U(0, 1) regardless of scale; reference behaviour
};

struct Dynamics{
    size_t pop_size = 20; // Constraint: > 0
    Real c1 = 2.0; // Constraint: >= 0
    Real c2 = 1.7; // Constraint: >= 0
    Real w_max = 0.5; // Constraint: > 0, > w_min
    Real w_min = 0.01; // Constraint: >= 0
    Real w_decay_frac = 0.75; // bounded in (0, 1]
    // Fraction of max_iter over which w decays w_max -> w_min, then holds.
    InertiaSchedule w_schedule = InertiaSchedule::Linear;

    Real vclamp_frac = 0.5; // Constraint: > 0
    // Per-dimension speed cap as a fraction of that dimension's span:
    // v_max[i] = vclamp_frac * (upper[i] - lower[i]).

    VelocityInit velocity_init = VelocityInit::ScaledUniform;

    bool constriction = false;
    // Clerc factor chi = |2k / (2 - phi - sqrt(phi^2 - 4 phi))|, phi = c1 + c2,
    // applied after the velocity update. Constraint if enabled: c1 + c2 > 4.
    // Double-damping: enabling this alongside a decaying w is usually a mistake.

    Real constriction_k = 1.0; // Constraint if enabled: != 0
};

// ---------------------------------------------------------------------------
// Topology (source of the social attractor)
// ---------------------------------------------------------------------------

// Which particles a given particle may learn from, i.e. where the c2 social
// attractor comes from. Global converges fastest but collapses the swarm onto
// one basin; the neighbourhood kinds trade convergence speed for the ability
// to hold several basins at once.
enum class TopologyKind : uint8_t{
    Global, // single swarm best
    Ring,   // static index-adjacent neighbourhood
    Knn,    // exact k nearest by L2, O(pop_size^2)
    Lsh     // approximate k nearest by locality-sensitive hashing
};

struct LshParams{
    size_t n_projections = 5; // k, constraint: >= 4
    // Number of h() projections composing one hash g(). Higher = more selective.
    Real window = 3.0; // w, constraint: >= 2
    // Quantisation window of h(x) = floor((v.x + b) / w). Higher = coarser.
    size_t n_tables = 30; // L, constraint: >= 5
    // Independent hash tables. Higher = better recall, cost is linear in L.
    size_t n_buckets = 0; // 0 = auto (pop_size / 2, minimum 1)
};

struct Topology{
    TopologyKind kind = TopologyKind::Global;
    size_t n_neighbours = 0; // 0 = auto (pop_size / 4, minimum 1)
    // Constraint: n_neighbours <= pop_size.
    Real radius_frac = 1.0; // bounded in (0, 1]
    // Candidates restricted to a hypercube of side radius_frac * (upper - lower)
    // centred on the query particle. 1.0 = no spatial restriction.
    Real dedup_tol = 1e-12; // Constraint: >= 0
    // Two particles closer than this count as the same point when collecting
    // neighbours. Distinct from SwarmConvergence::objfun_tol.
    LshParams lsh;

    // The social attractor is selected on personal-best VALUES, and falls back
    // to the swarm best when a neighbourhood is empty or entirely infinite.
};


// ---------------------------------------------------------------------------
// Multi-minimum strategy
// ---------------------------------------------------------------------------

// The outer loop wrapping the swarm: how many swarms run, over what domains,
// and what makes them converge to DIFFERENT minima. All three reduce to
// repeated invocations of the same core swarm.
enum class SearchStrategy : uint8_t{
    Single,     // one swarm, one minimum
    Deflection, // sequential restarts, each repelled from prior minima
    Grid        // independent swarms over a subdivision of the domain
};

struct Deflection{
    size_t runs = 1; // Constraint: >= 1
    Real lambda = 50.0; // Constraint: > 0
    // Sharpness of the pole: T(x) = tanh(lambda * ||x - x*||). Large lambda
    // gives a narrow spike; small lambda a wide repulsive region.
    Real shift = 0.0;
    // F(x) = (f(x) + shift) * prod_i 1 / T_i(x). Keeps F positive when the
    // target value is not 0, so the poles actually repel.
    bool stop_when_barren = true;
    // End early once a run yields no new minimum.
};

struct Repulsion{
    bool enabled = false;
    Real radius = 0.25; // Constraint if enabled: > 0
    Real rho = 3.0; // Constraint if enabled: != 0
    // Geometric kick along the outward unit vector, applied after the position
    // update. Independent of and composable with Deflection, which is functional.
};

struct Grid{
    size_t resolution = 1; // Constraint: >= 0
    // Each dimension is split into 2^resolution slices, giving
    // 2^(PSO_DIM * resolution) subspaces.
    size_t start_subspace = 0; // Constraint: < 2^(PSO_DIM * resolution)
    size_t batch = 0; // 0 = auto (hardware threads / execution.threads)
    // Concurrent subspace workers, drawn from a fixed pool.
};

// ---------------------------------------------------------------------------
// Objective
// ---------------------------------------------------------------------------

// Which objective plugin to evaluate. Selects both f(x) and which parameter
// sub-struct below is read; every other block is objective-agnostic.
enum class ObjectiveKind : uint8_t{
    Poincare,   // distance from a periodic orbit of a caldera Hamiltonian
    Analytic    // standard test functions, for exercising the swarm itself
};

// How "the trajectory returned to where it started" is reduced to the single
// number the swarm minimises, once the return point has been integrated.
enum class PoincareMetric : uint8_t{
    Weighted,  // magnitude-weighted, normalised residual
    Euclidean  // plain L2 return distance
};

struct PoincareParams{
    Real energy = 17.0;
    // Conserved energy E. Defines the accessible region; particles off the
    // energy surface evaluate to +inf.
    Real section = -1.8019693;
    // The y value defining the surface of section {y = section}. A particle
    // (x, px) is lifted to q = (x, section, px, py), py solved from E.
    size_t n_crossings = 1; // Constraint: >= 1
    // Crossings before comparing back to the start. 1 = fixed point of the
    // return map, 2 = period-2 orbit. Also sets how many columns one recorded
    // minimum occupies.

    Real c1 = 5.0;    // quadratic term    c1 (x^2 + y^2)
    Real c2 = 3.0;    // linear tilt       c2 y
    Real c3 = -0.3;   // quartic term     -c3 (x^4 + y^4 - 6 x^2 y^2)

    std::string integrator = "RK45";
    // Stepping method, mapped onto ode::Stepper at solver construction.
    // Valid: "Euler", "RK4", "RK23", "RK45", "DOP853", "BDF". Validated at
    // startup, since an unrecognised name cannot be caught by the compiler.
    // Euler and RK4 are fixed-step and ignore rtol/atol. RK45 is the default
    // adaptive choice; DOP853 is higher order and takes fewer steps at the very
    // tight tolerances the LongReal phase needs. BDF targets stiff systems, so
    // it is unlikely to help for a Hamiltonian flow.

    Real dt_init = 1e-3; // Constraint: > 0
    // Initial step size. For adaptive methods only a starting guess; for the
    // fixed-step methods it is the step size for the whole integration.
    Real rtol = 1e-12; // Constraint: > 0
    // Integrator relative tolerance. Constraint: rtol <= 1e-3 * objfun_tol.
    // The objective cannot be more accurate than the integrator computing it.
    Real atol = 0.0; // Constraint: >= 0
    Real section_ftol = 0.0; // Constraint: >= 0
    // Root-find tolerance on the crossing condition y - section = 0.
    Real energy_drift_tol = 1e-6; // Constraint: > 0
    // Abort an integration once |E(t) - E| / E exceeds this.
    Real t_max = 1000.0; // Constraint: > 0

    PoincareMetric metric = PoincareMetric::Weighted;
};

// Parameters of the analytic test functions.
struct AnalyticParams{
    std::string function = "rastrigin";
    // Which function to minimise: "sphere", "rastrigin" or "rosenbrock".
    // All three have their global minimum at value 0, so objfun_target = 0
    // suits every one of them.
};

struct Objective{
    ObjectiveKind kind = ObjectiveKind::Poincare;
    // Each objective owns its own block below; only the selected one is read.
    PoincareParams poincare;
    AnalyticParams analytic;
};

// ---------------------------------------------------------------------------
// Execution
// ---------------------------------------------------------------------------

// Verbosity of per-run diagnostic output.
enum class LogLevel : uint8_t{
    Quiet,
    Info,
    Trace
};

struct Execution{
    uint64_t seed = 42;
    // Master seed. Each worker draws from mix(seed, worker_id) so that runs stay
    // reproducible while subspaces explore independently.
    size_t threads = 1; // Constraint: >= 1
    // Objective-evaluation parallelism; particles are embarrassingly parallel.
    std::string output_dir = "./out";
    bool resume = false;
    // When false the output directory is cleared on startup. Must be true to
    // continue a run via Grid::start_subspace.
    LogLevel log_level = LogLevel::Info;
};

// ---------------------------------------------------------------------------
// Aggregate
// ---------------------------------------------------------------------------

struct Parameters{
    Domain domain;
    SwarmConvergence convergence;
    Dynamics dynamics;
    Topology topology;

    SearchStrategy strategy = SearchStrategy::Single;
    Deflection deflection;
    Repulsion repulsion;
    Grid grid;

    Objective objective;
    Execution execution;
};

} // namespace input

#endif // PSO_INPUT_PARAMETERS_HPP
