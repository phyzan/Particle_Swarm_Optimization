#ifndef PSO_SWARM_TOPOLOGY_HPP
#define PSO_SWARM_TOPOLOGY_HPP

#include <src/input/parameters.hpp>
#include <algorithm>
#include <array>
#include <random>
#include <unordered_map>
#include <vector>

using std::abs;

namespace topology{

// Candidate neighbours of one particle, by index. Whoever calls this picks
// the best of them by personal-best VALUE, not by distance: the distance only
// decides who is in the neighbourhood.
//
// Every kind must tolerate returning nothing. An empty result means "no
// neighbourhood", and the caller falls back to the swarm best rather than
// reading an uninitialised attractor.
class Neighbourhood{
public:

    /// @brief Builds the neighbourhood policy the parameters select.
    Neighbourhood(const input::Topology& tp, size_t pop_size,
                  const std::array<Real, PSO_DIM>& span);

    /// @brief Which neighbourhood kind is in use.
    input::TopologyKind mode() const;

    /// @brief Rebuilds the per-iteration hash tables; a no-op unless kind is Lsh.
    template<typename Positions>
    void rebuild(const Positions& X, size_t pop_size, std::mt19937_64& rng)
    {
        if (kind != input::TopologyKind::Lsh){
            return;
        }

        projections.assign(n_tables*n_projections, std::array<Real, PSO_DIM>{});
        offsets.assign(n_tables*n_projections, 0.0);
        mixers.assign(n_tables*n_projections, 0);

        std::normal_distribution<Real> gauss(0.0, 1.0);
        std::uniform_real_distribution<Real> offset(0.0, window);
        std::uniform_int_distribution<uint64_t> mix(1, MOD - 1);

        for (size_t s = 0; s < n_tables*n_projections; s++){
            for (size_t i = 0; i < PSO_DIM; i++){
                projections[s][i] = gauss(rng);
            }
            offsets[s] = offset(rng);
            mixers[s] = mix(rng);
        }

        tables.assign(n_tables, {});

        for (size_t t = 0; t < n_tables; t++){
            for (size_t j = 0; j < pop_size; j++){
                tables[t].emplace(slot(t, X, j), j);
            }
        }
    }

    /// @brief Collects the neighbour indices of particle j into `out`.
    template<typename Positions>
    void query(const Positions& X, size_t j, size_t pop_size, std::vector<size_t>& out) const
    {
        out.clear();

        switch (kind){

            case input::TopologyKind::Global:
                return;

            case input::TopologyKind::Ring:
                out.push_back((j + pop_size - 1) % pop_size);
                out.push_back(j);
                out.push_back((j + 1) % pop_size);
                return;

            case input::TopologyKind::Knn:
                for (size_t q = 0; q < pop_size; q++){
                    if (q != j && within_box(X, j, q)){
                        out.push_back(q);
                    }
                }
                break;

            case input::TopologyKind::Lsh:
                for (size_t t = 0; t < n_tables; t++){
                    const auto range = tables[t].equal_range(slot(t, X, j));
                    for (auto it = range.first; it != range.second; ++it){
                        if (it->second != j && within_box(X, j, it->second)){
                            out.push_back(it->second);
                        }
                    }
                }
                std::sort(out.begin(), out.end());
                out.erase(std::unique(out.begin(), out.end()), out.end());
                break;
        }

        drop_duplicates(X, j, out);
        keep_nearest(X, j, out);
    }

private:

    static constexpr uint64_t MOD = (uint64_t{1} << 61) - 1;

    /// @brief Whether q lies inside j's neighbourhood box.
    template<typename Positions>
    bool within_box(const Positions& X, size_t j, size_t q) const{
        for (size_t i = 0; i < PSO_DIM; i++){
            if (abs(Real(X(q, i)) - Real(X(j, i))) > half_box[i]){
                return false;
            }
        }
        return true;
    }

    /// @brief Euclidean distance between particles j and q.
    template<typename Positions>
    Real distance(const Positions& X, size_t j, size_t q) const{
        Real acc = 0;
        for (size_t i = 0; i < PSO_DIM; i++){
            const Real d = Real(X(q, i)) - Real(X(j, i));
            acc += d*d;
        }
        return std::sqrt(acc);
    }

    /// @brief Removes candidates closer to j than dedup_tol.
    template<typename Positions>
    void drop_duplicates(const Positions& X, size_t j, std::vector<size_t>& out) const{
        out.erase(std::remove_if(out.begin(), out.end(),
                                 [&](size_t q){ return distance(X, j, q) <= dedup_tol; }),
                  out.end());
    }

    /// @brief Trims the candidates to the k nearest to j.
    template<typename Positions>
    void keep_nearest(const Positions& X, size_t j, std::vector<size_t>& out) const{
        if (out.size() <= k){
            return;
        }
        std::nth_element(out.begin(), out.begin() + std::ptrdiff_t(k), out.end(),
                         [&](size_t a, size_t b){ return distance(X, j, a) < distance(X, j, b); });
        out.resize(k);
    }

    /// @brief Hashes particle j into a bucket of the given table.
    template<typename Positions>
    uint64_t slot(size_t table, const Positions& X, size_t j) const{

        // h(x) = floor((a.x + b)/w), kept SIGNED: a.x is Gaussian and
        // routinely negative, and an unsigned cast here destroys the hash.
        // Each term is reduced modulo MOD before accumulating so nothing
        // overflows.
        uint64_t id = 0;

        for (size_t m = 0; m < n_projections; m++){

            const size_t s = table*n_projections + m;

            Real dot = offsets[s];
            for (size_t i = 0; i < PSO_DIM; i++){
                dot += projections[s][i] * Real(X(j, i));
            }

            const auto h = int64_t(std::floor(dot / window));
            const auto folded = uint64_t(h % int64_t(MOD) + int64_t(MOD)) % MOD;

            id = (id + (mixers[s] % MOD) * folded) % MOD;
        }

        return id % n_buckets;
    }

    input::TopologyKind kind;
    size_t k;
    Real dedup_tol;
    std::array<Real, PSO_DIM> half_box{};

    size_t n_tables;
    size_t n_projections;
    Real window;
    size_t n_buckets;

    std::vector<std::array<Real, PSO_DIM>> projections;
    std::vector<Real> offsets;
    std::vector<uint64_t> mixers;
    std::vector<std::unordered_multimap<uint64_t, size_t>> tables;
};

} // namespace topology

#endif // PSO_SWARM_TOPOLOGY_HPP
