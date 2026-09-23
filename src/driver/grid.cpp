#include <src/driver/grid.hpp>

#include <src/driver/search.hpp>
#include <src/objective/objective.hpp>

#include <algorithm>
#include <atomic>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <thread>
#include <vector>

namespace driver{

/// @brief Number of subspaces the resolution produces.
size_t grid_size(const input::Grid& cfg)
{
    return size_t{1} << (PSO_DIM * cfg.resolution);
}

/// @brief Bounds of one subspace, by mixed-radix decomposition of its index.
void subspace_bounds(size_t cell, const input::Parameters& params,
                            std::array<Real, PSO_DIM>& lower,
                            std::array<Real, PSO_DIM>& upper)
{
    // Computed on demand: a resolution-4 grid in 2D is 256 cells, and there
    // is no reason to materialise them.
    const size_t divisions = size_t{1} << params.grid.resolution;

    for (size_t i = 0; i < PSO_DIM; i++){

        const Real step = (params.domain.upper[i] - params.domain.lower[i]) / Real(divisions);

        size_t stride = 1;
        for (size_t k = i + 1; k < PSO_DIM; k++){
            stride *= divisions;
        }

        const size_t digit = (cell / stride) % divisions;

        lower[i] = params.domain.lower[i] + Real(digit)*step;
        upper[i] = lower[i] + step;
    }
}

/// @brief Derives an independent RNG stream for one worker.
uint64_t seed_for(uint64_t master, uint64_t worker)
{
    auto mix = [](uint64_t z){
        z += 0x9E3779B97F4A7C15ULL;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    };
    return mix(master ^ mix(worker));
}


Minima run_grid(const input::Parameters& params, std::ostream& out)
{
    namespace fs = std::filesystem;

    const size_t cells = grid_size(params.grid);
    const fs::path destination(params.execution.output_dir);

    // Only clear when this is not a resume: validate() already refuses
    // start_subspace > 0 without resume, which is how the reference
    // implementation destroyed the results it was meant to continue.
    if (!params.execution.resume && fs::exists(destination)){
        fs::remove_all(destination);
    }

    fs::create_directories(destination);

    // Subspaces are independent and evenly loaded, so batch parallelises far
    // better than execution.threads does -- that one is capped by a 20-particle
    // loop with wildly uneven work. But the two MULTIPLY, so auto splits the
    // machine between them rather than letting them oversubscribe it.
    const size_t hardware = std::max<size_t>(1, std::thread::hardware_concurrency());
    const size_t requested = params.grid.batch
                           ? params.grid.batch
                           : std::max<size_t>(1, hardware / std::max<size_t>(1, params.execution.threads));

    const size_t workers = std::min<size_t>(requested, cells - params.grid.start_subspace);

    out << "grid: " << cells << " subspaces, " << workers << " at a time";

    if (params.grid.start_subspace > 0){
        out << ", resuming from " << params.grid.start_subspace;
    }

    out << "\n";

    Minima found;
    std::mutex guard;
    std::atomic<size_t> next{params.grid.start_subspace};

    const auto worker = [&]{

        for (size_t cell = next++; cell < cells; cell = next++){

            input::Parameters local = params;

            subspace_bounds(cell, params, local.domain.lower, local.domain.upper);
            local.execution.seed = seed_for(params.execution.seed, cell);

            std::ostringstream name;
            name << "subspace_" << std::setw(5) << std::setfill('0') << cell << ".txt";

            std::ofstream file(destination / name.str());
            file << std::setprecision(16);

            file << "subspace " << cell << " of " << cells << '\n'
                 << "bounds  [" << local.domain.lower[0] << ", " << local.domain.upper[0] << "]"
                 << " x [" << local.domain.lower[1] << ", " << local.domain.upper[1] << "]\n";

            Minima local_found = run_swarm(local, file);

            {
                std::lock_guard<std::mutex> lock(guard);

                found.insert(found.end(), local_found.begin(), local_found.end());

                out << "  subspace " << std::setw(5) << std::setfill('0') << cell
                    << " done (" << local_found.size()/objective_stride(params)
                    << " minima)\n";
            }
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(workers);

    for (size_t w = 0; w < workers; w++){
        pool.emplace_back(worker);
    }

    for (auto& thread : pool){
        thread.join();
    }

    return found;
}

} // namespace driver
