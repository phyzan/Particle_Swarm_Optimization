#include <src/input/initialize.hpp>

namespace input{

/// @brief Parses a TOML file from disk.
InputFile::InputFile(const std::string& path) : tbl_(toml::parse_file(path)) {}

/// @brief Wraps an already-parsed TOML table.
InputFile::InputFile(toml::table table) : tbl_(std::move(table)) {}


/// @brief Reads a high-precision scalar, keeping every digit of a quoted literal.
LongReal InputFile::get_or_add_long(std::string_view block, std::string_view key, LongReal fallback)
{
    const auto node = mark(block, key);

    if (!node){
        return fallback;
    }

    if (const auto text = node.template value<std::string>()){
#if PSO_USE_MPREAL_SWAP
        return LongReal(text->c_str(), mpfr::mpreal::get_default_prec());
#else
        return LongReal(std::stold(*text));
#endif
    }

    const auto value = node.template value<Real>();

    if (!value){
        throw std::runtime_error(location(node) + ": '" + name(block, key)
                               + "' must be a number or a quoted number");
    }

    return LongReal(*value);
}


/// @brief Lists file keys no accessor asked for, with their line numbers.
std::vector<std::string> InputFile::unknown_keys() const
{
    std::vector<std::string> unknown;

    for (const auto& [block, node] : tbl_){
        const auto* table = node.as_table();

        if (!table){
            unknown.emplace_back(std::string(block.str())
                               + "  (line " + std::to_string(node.source().begin.line)
                               + ", not inside a [block])");
            continue;
        }

        for (const auto& [key, value] : *table){
            const std::string full = std::string(block.str()) + '.' + std::string(key.str());

            if (!touched_.contains(full)){
                unknown.emplace_back(full + "  (line "
                                   + std::to_string(value.source().begin.line) + ")");
            }
        }
    }

    return unknown;
}


/// @brief Binds an input file onto Parameters, field by field.
Parameters read_parameters(InputFile& in)
{
// Each line passes the struct's own member as the fallback, so defaults
// are declared once in parameters.hpp and cannot drift away from it.
Parameters p{};

// ---- [domain] -- lower/upper have no defaults, so they are required ----
p.domain.lower = in.get_array<PSO_DIM>("domain", "lower");
p.domain.upper = in.get_array<PSO_DIM>("domain", "upper");
p.domain.bound_policy = in.get_or_add_enum("domain", "bound_policy", p.domain.bound_policy,
    {{"none", BoundPolicy::None}, {"clamp", BoundPolicy::Clamp},
     {"reflect", BoundPolicy::Reflect}, {"reinit", BoundPolicy::ReInit}});

// ---- [convergence] ----
auto& cv = p.convergence;
cv.objfun_target   = in.get_or_add<Real>  ("convergence", "objfun_target",   cv.objfun_target);
cv.objfun_tol      = in.get_or_add<Real>  ("convergence", "objfun_tol",      cv.objfun_tol);
cv.swap_tol        = in.get_or_add<Real>  ("convergence", "swap_tol",        cv.swap_tol);
cv.mpfr_prec       = in.get_or_add<int>   ("convergence", "mpfr_prec",       cv.mpfr_prec);
cv.max_iter        = in.get_or_add<size_t>("convergence", "max_iter",        cv.max_iter);
cv.max_stall_iter  = in.get_or_add<size_t>("convergence", "max_stall_iter",  cv.max_stall_iter);

// ---- [dynamics] ----
auto& dy = p.dynamics;
dy.pop_size       = in.get_or_add<size_t>("dynamics", "pop_size",       dy.pop_size);
dy.c1             = in.get_or_add<Real>  ("dynamics", "c1",             dy.c1);
dy.c2             = in.get_or_add<Real>  ("dynamics", "c2",             dy.c2);
dy.w_max          = in.get_or_add<Real>  ("dynamics", "w_max",          dy.w_max);
dy.w_min          = in.get_or_add<Real>  ("dynamics", "w_min",          dy.w_min);
dy.w_decay_frac   = in.get_or_add<Real>  ("dynamics", "w_decay_frac",   dy.w_decay_frac);
dy.vclamp_frac    = in.get_or_add<Real>  ("dynamics", "vclamp_frac",    dy.vclamp_frac);
dy.constriction   = in.get_or_add<bool>  ("dynamics", "constriction",   dy.constriction);
dy.constriction_k = in.get_or_add<Real>  ("dynamics", "constriction_k", dy.constriction_k);
dy.w_schedule = in.get_or_add_enum("dynamics", "w_schedule", dy.w_schedule,
    {{"linear", InertiaSchedule::Linear}, {"exponential", InertiaSchedule::Exponential},
     {"constant", InertiaSchedule::Constant}});
dy.velocity_init = in.get_or_add_enum("dynamics", "velocity_init", dy.velocity_init,
    {{"zero", VelocityInit::Zero}, {"scaled_uniform", VelocityInit::ScaledUniform},
     {"uniform01", VelocityInit::Uniform01}});

// ---- [topology] ----
auto& tp = p.topology;
tp.n_neighbours = in.get_or_add<size_t>("topology", "n_neighbours", tp.n_neighbours);
tp.radius_frac  = in.get_or_add<Real>  ("topology", "radius_frac",  tp.radius_frac);
tp.dedup_tol    = in.get_or_add<Real>  ("topology", "dedup_tol",    tp.dedup_tol);
tp.kind = in.get_or_add_enum("topology", "kind", tp.kind,
    {{"global", TopologyKind::Global}, {"ring", TopologyKind::Ring},
     {"knn", TopologyKind::Knn}, {"lsh", TopologyKind::Lsh}});

// ---- [lsh] ----
auto& ls = p.topology.lsh;
ls.n_projections = in.get_or_add<size_t>("lsh", "n_projections", ls.n_projections);
ls.window        = in.get_or_add<Real>  ("lsh", "window",        ls.window);
ls.n_tables      = in.get_or_add<size_t>("lsh", "n_tables",      ls.n_tables);
ls.n_buckets     = in.get_or_add<size_t>("lsh", "n_buckets",     ls.n_buckets);

// ---- [strategy] ----
p.strategy = in.get_or_add_enum("strategy", "kind", p.strategy,
    {{"single", SearchStrategy::Single}, {"deflection", SearchStrategy::Deflection},
     {"grid", SearchStrategy::Grid}});

// ---- [deflection] ----
auto& df = p.deflection;
df.runs              = in.get_or_add<size_t>("deflection", "runs",              df.runs);
df.lambda            = in.get_or_add<Real>  ("deflection", "lambda",            df.lambda);
df.shift             = in.get_or_add<Real>  ("deflection", "shift",             df.shift);
df.stop_when_barren  = in.get_or_add<bool>  ("deflection", "stop_when_barren",  df.stop_when_barren);

// ---- [repulsion] ----
auto& rp = p.repulsion;
rp.enabled = in.get_or_add<bool>("repulsion", "enabled", rp.enabled);
rp.radius  = in.get_or_add<Real>("repulsion", "radius",  rp.radius);
rp.rho     = in.get_or_add<Real>("repulsion", "rho",     rp.rho);

// ---- [grid] ----
auto& gr = p.grid;
gr.resolution      = in.get_or_add<size_t>("grid", "resolution",      gr.resolution);
gr.start_subspace  = in.get_or_add<size_t>("grid", "start_subspace",  gr.start_subspace);
gr.batch           = in.get_or_add<size_t>("grid", "batch",           gr.batch);

// ---- [objective] ----
p.objective.kind = in.get_or_add_enum("objective", "kind", p.objective.kind,
    {{"poincare", ObjectiveKind::Poincare}, {"analytic", ObjectiveKind::Analytic}});

// ---- [poincare] ----
auto& pc = p.objective.poincare;
pc.energy           = in.get_or_add<Real>       ("poincare", "energy",           pc.energy);
pc.section          = in.get_or_add<Real>       ("poincare", "section",          pc.section);
pc.n_crossings      = in.get_or_add<size_t>     ("poincare", "n_crossings",      pc.n_crossings);
pc.c1               = in.get_or_add<Real>       ("poincare", "c1",               pc.c1);
pc.c2               = in.get_or_add<Real>       ("poincare", "c2",               pc.c2);
pc.c3               = in.get_or_add<Real>       ("poincare", "c3",               pc.c3);
pc.integrator       = in.get_or_add<std::string>("poincare", "integrator",       pc.integrator);
pc.dt_init          = in.get_or_add<Real>       ("poincare", "dt_init",          pc.dt_init);
pc.rtol             = in.get_or_add<Real>       ("poincare", "rtol",             pc.rtol);
pc.atol             = in.get_or_add<Real>       ("poincare", "atol",             pc.atol);
pc.section_ftol     = in.get_or_add<Real>       ("poincare", "section_ftol",     pc.section_ftol);
pc.energy_drift_tol = in.get_or_add<Real>       ("poincare", "energy_drift_tol", pc.energy_drift_tol);
pc.t_max            = in.get_or_add<Real>       ("poincare", "t_max",            pc.t_max);
pc.metric = in.get_or_add_enum("poincare", "metric", pc.metric,
    {{"weighted", PoincareMetric::Weighted}, {"euclidean", PoincareMetric::Euclidean}});

// ---- [analytic] ----
p.objective.analytic.function =
    in.get_or_add<std::string>("analytic", "function", p.objective.analytic.function);

// ---- [execution] ----
auto& ex = p.execution;
ex.seed       = in.get_or_add<uint64_t>   ("execution", "seed",       ex.seed);
ex.threads    = in.get_or_add<size_t>     ("execution", "threads",    ex.threads);
ex.output_dir = in.get_or_add<std::string>("execution", "output_dir", ex.output_dir);
ex.resume     = in.get_or_add<bool>       ("execution", "resume",     ex.resume);
ex.log_level  = in.get_or_add_enum("execution", "log_level", ex.log_level,
    {{"quiet", LogLevel::Quiet}, {"info", LogLevel::Info}, {"trace", LogLevel::Trace}});

return p;
}


/// @brief Reads a file and rejects it if it names anything unknown.
Parameters read_parameters(const std::string& path)
{
InputFile in(path);

Parameters p = read_parameters(in);

const auto unknown = in.unknown_keys();

if (!unknown.empty()){
    std::string message = "unrecognised parameters in " + path + ":";

    for (const auto& key : unknown){
        message += "\n  " + key;
    }

    throw std::runtime_error(message);
}

return p;
}


} // namespace input
