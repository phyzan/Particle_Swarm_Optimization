#include <src/swarm/topology.hpp>

namespace topology{

/// @brief Builds the neighbourhood policy the parameters select.
Neighbourhood::Neighbourhood(const input::Topology& tp, size_t pop_size,
                       const std::array<Real, PSO_DIM>& span)
    : kind(tp.kind),
      k(tp.n_neighbours ? tp.n_neighbours : std::max<size_t>(1, pop_size/4)),
      dedup_tol(tp.dedup_tol),
      n_tables(tp.lsh.n_tables),
      n_projections(tp.lsh.n_projections),
      window(tp.lsh.window),
      n_buckets(tp.lsh.n_buckets ? tp.lsh.n_buckets : std::max<size_t>(1, pop_size/2))
{
    for (size_t i = 0; i < PSO_DIM; i++){
        half_box[i] = tp.radius_frac * span[i] / 2;
    }

    k = std::min(k, pop_size);
}

/// @brief Which neighbourhood kind is in use.
input::TopologyKind Neighbourhood::mode() const
{
    return kind;
}


} // namespace topology
