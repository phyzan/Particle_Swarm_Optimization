#ifndef PSO_SWARM_TOPOLOGY_HPP
#define PSO_SWARM_TOPOLOGY_HPP

// #include "swarm/swarm.hpp"
#include <memory>
#include <src/input/parameters.hpp>
#include <algorithm>
#include <array>
#include <random>
#include <unordered_map>
#include <vector>
#include <ndspan/ndspan.hpp>

using std::abs;
using ndspan::Array2D;

namespace topology{

// Candidate neighbours of one particle, by index. Whoever calls this picks
// the best of them by personal-best VALUE, not by distance: the distance only
// decides who is in the neighbourhood.
//
// Every kind must tolerate returning nothing. An empty result means "no
// neighbourhood", and the caller falls back to the swarm best rather than
// reading an uninitialised attractor.
//
// Index convention used throughout, in signatures and comments alike: j is the
// query particle, the one whose neighbourhood is being built, and q is a
// candidate neighbour of it.
template<typename T>
class Neighbourhood{
public:

    /// @brief Builds the neighbourhood policy the parameters select.
    Neighbourhood(const input::Topology& tp, size_t pop_size,
                        const std::array<Real, PSO_DIM>& span)
        : pop_size(pop_size),
        k(tp.n_neighbours ? tp.n_neighbours : std::max<size_t>(1, pop_size/4)),
        dedup_tol(tp.dedup_tol){
        for (size_t i = 0; i < PSO_DIM; i++){
            half_box[i] = tp.radius_frac * span[i] / 2;
        }

        k = std::min(k, pop_size);
    }

    /// @brief Rebuilds the per-iteration hash tables; a no-op by default
    virtual void rebuild(const Array2D<T, 0, PSO_DIM>& X, std::mt19937_64& rng) {}

    /// @brief Collects the neighbour indices of particle j into `out`.
    virtual std::vector<size_t> query(const Array2D<T, 0, PSO_DIM>& X, size_t j) const {
        return {};
    }

protected:

    /// @brief Whether q lies inside j's neighbourhood box.
    bool within_box(const Array2D<T, 0, PSO_DIM>& X, size_t j, size_t q) const{
        for (size_t i = 0; i < PSO_DIM; i++){
            if (abs(X(q, i) - X(j, i)) > half_box[i]){
                return false;
            }
        }
        return true;
    }

    /// @brief Euclidean distance between particles j and q.
    T distance(const Array2D<T, 0, PSO_DIM>& X, size_t j, size_t q) const{
        T acc = 0;
        for (size_t i = 0; i < PSO_DIM; i++){
            const auto d = X(q, i) - X(j, i);
            acc += d*d;
        }
        return sqrt(acc);
    }

    /// @brief Removes candidates closer to j than dedup_tol.
    void drop_duplicates(const Array2D<T, 0, PSO_DIM>& X, size_t j, std::vector<size_t>& out) const{
        out.erase(std::remove_if(out.begin(), out.end(),
                                 [&](size_t q){ return distance(X, j, q) <= dedup_tol; }),
                  out.end());
    }

    /// @brief Trims the candidates to the k nearest to j.
    void keep_nearest(const Array2D<T, 0, PSO_DIM>& X, size_t j, std::vector<size_t>& out) const{
        if (out.size() <= k){
            return;
        }
        std::nth_element(out.begin(), out.begin() + std::ptrdiff_t(k), out.end(),
                         [&](size_t a, size_t b){ return distance(X, j, a) < distance(X, j, b); });
        out.resize(k);
    }

    size_t pop_size;
    size_t k;                   // most neighbours kept, nearest first (0 in => pop_size/4)
    Real dedup_tol;             // a candidate this close to the query particle is the same point, so it goes
    std::array<T, PSO_DIM> half_box{};
    // Half-width of j's candidate box, per dimension: radius_frac * span[i] / 2. within_box()
    // tests every dimension against it, so the neighbourhood is a hypercube, not a ball.

};

template<typename T, input::TopologyKind TK>
class DerivedNeighbourhood;

template<typename T>
class DerivedNeighbourhood<T, input::TopologyKind::Lsh> : public Neighbourhood<T>{
    using Base = Neighbourhood<T>;
    static constexpr uint64_t MOD = (uint64_t{1} << 61) - 1;
public:
    /// @brief Builds the neighbourhood policy the parameters select.
    DerivedNeighbourhood(const input::Topology& tp, size_t pop_size,
                        const std::array<Real, PSO_DIM>& span)
        : Base(tp, pop_size, span),
        n_tables(tp.lsh.n_tables),
        n_projections(tp.lsh.n_projections),
        window(tp.lsh.window),
        n_buckets(tp.lsh.n_buckets ? tp.lsh.n_buckets : std::max<size_t>(1, pop_size/2)) {}

    /// @brief Rebuilds the per-iteration hash tables; a no-op by default
    void rebuild(const Array2D<T, 0, PSO_DIM>& X, std::mt19937_64& rng) override{
        projections.assign(n_tables*n_projections, std::array<T, PSO_DIM>{});
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
            for (size_t j = 0; j < this->pop_size; j++){
                tables[t].emplace(slot(t, X, j), j);
            }
        }
    }

    /// @brief Collects the neighbour indices of particle j into `out`.
    std::vector<size_t> query(const Array2D<T, 0, PSO_DIM>& X, size_t j) const override {
        std::vector<size_t> out;
        
        for (size_t t = 0; t < n_tables; t++){
            const auto range = tables[t].equal_range(slot(t, X, j));
            for (auto it = range.first; it != range.second; ++it){
                if (it->second != j && this->within_box(X, j, it->second)){
                    out.push_back(it->second);
                }
            }
        }
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());

        this->drop_duplicates(X, j, out);
        this->keep_nearest(X, j, out);

        return out;
    }

private:

    /// @brief Hashes particle j into a bucket of the given table.
    uint64_t slot(size_t table, const Array2D<T, 0, PSO_DIM>& X, size_t j) const{

        // h(x) = floor((a.x + b)/w), kept SIGNED: a.x is Gaussian and
        // routinely negative, and an unsigned cast here destroys the hash.
        // Each term is reduced modulo MOD before accumulating so nothing
        // overflows.
        uint64_t id = 0;

        for (size_t m = 0; m < n_projections; m++){

            const size_t s = table*n_projections + m;

            T dot = offsets[s];
            for (size_t i = 0; i < PSO_DIM; i++){
                dot += projections[s][i] * X(j, i);
            }

            const auto h = int64_t(std::floor(Real(dot) / window));
            const auto folded = uint64_t(h % int64_t(MOD) + int64_t(MOD)) % MOD;

            id = (id + (mixers[s] % MOD) * folded) % MOD;
        }

        return id % n_buckets;
    }
    // ---- LSH hyperparameters. Read only when kind == Lsh; the usual L / k / w of p-stable LSH.
    size_t n_tables;            // L: independent tables, whose buckets query() unions. Recall vs cost, linear
    size_t n_projections;       // k: h() projections ANDed into one table's hash. Higher = more selective
    Real window;                // w: quantisation width of h(). Higher = coarser buckets, more collisions
    size_t n_buckets;           // buckets per table, the modulus slot() ends on (0 in => pop_size/2)

    // ---- LSH state, discarded and redrawn by every rebuild() since the particles have moved.
    // All three are flat [n_tables * n_projections], indexed s = table*n_projections + m.
    std::vector<std::array<T, PSO_DIM>> projections;  // a: N(0,1) per component -- Gaussian is what
                                                         // makes h() distance-preserving in L2
    std::vector<T> offsets;                           // b: uniform [0, window), so bucket edges do not
                                                         // sit at fixed positions
    std::vector<uint64_t> mixers;                        // random in [1, MOD), folding the n_projections
                                                         // values into one id. Hash mixing, not LSH
    std::vector<std::unordered_multimap<uint64_t, size_t>> tables;  // per table: bucket id -> particle index
};


template<typename T>
class DerivedNeighbourhood<T, input::TopologyKind::Ring> : public Neighbourhood<T>{
    using Base = Neighbourhood<T>;
    
public:
    using Base::Base;

    /// @brief Collects the neighbour indices of particle j into `out`.
    std::vector<size_t> query(const Array2D<T, 0, PSO_DIM>& X, size_t j) const override {
        std::vector<size_t> out{3};

        out[0] = (j + this->pop_size - 1) % this->pop_size;
        out[1] = j;
        out[2] = (j + 1) % this->pop_size;
        
        Base::drop_duplicates(X, j, out);
        Base::keep_nearest(X, j, out);

        return out;
    }
};


template<typename T>
class DerivedNeighbourhood<T, input::TopologyKind::Knn> : public Neighbourhood<T>{
    using Base = Neighbourhood<T>;
    
public:
    using Base::Base;

    /// @brief Collects the neighbour indices of particle j into `out`.
    std::vector<size_t> query(const Array2D<T, 0, PSO_DIM>& X, size_t j) const override {
        std::vector<size_t> out;

        for (size_t q = 0; q < this->pop_size; q++){
            if (q != j && Base::within_box(X, j, q)){
                out.push_back(q);
            }
        }

        Base::drop_duplicates(X, j, out);
        Base::keep_nearest(X, j, out);

        return out;
    }
};


template<typename T>
std::unique_ptr<Neighbourhood<T>> make_neighbourhood(const input::Topology& tp, size_t pop_size, const std::array<Real, PSO_DIM>& span){
    switch (tp.kind){
        case input::TopologyKind::Lsh:
            return std::make_unique<DerivedNeighbourhood<T, input::TopologyKind::Lsh>>(tp, pop_size, span);
        case input::TopologyKind::Ring:
            return std::make_unique<DerivedNeighbourhood<T, input::TopologyKind::Ring>>(tp, pop_size, span);
        case input::TopologyKind::Knn:
            return std::make_unique<DerivedNeighbourhood<T, input::TopologyKind::Knn>>(tp, pop_size, span);
        default:
            return std::make_unique<Neighbourhood<T>>(tp, pop_size, span);
    }
}


} // namespace topology

#endif // PSO_SWARM_TOPOLOGY_HPP
