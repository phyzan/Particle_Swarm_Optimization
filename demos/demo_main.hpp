#ifndef PSO_DEMOS_DEMO_MAIN_HPP
#define PSO_DEMOS_DEMO_MAIN_HPP

// Shared body for the demo binaries. Each demo is a label plus a .toml path
// baked in at configure time (PSO_DEMO_INPUT), so the binaries take no
// arguments and find their input from any working directory.
//
// Demos never hardcode parameters -- everything comes from the file.

#include <src/driver/dispatch.hpp>
#include <src/input/initialize.hpp>
#include <src/input/validate.hpp>
#include <src/objective/objective.hpp>

#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

namespace demo{

/// @brief Name of a search strategy, for the banner.
inline const char* strategy_name(input::SearchStrategy s)
{
    switch (s){
    case input::SearchStrategy::Single:     return "single";
    case input::SearchStrategy::Deflection: return "deflection";
    case input::SearchStrategy::Grid:       return "grid";
    }
    return "unknown";
}

/// @brief Name of a neighbourhood topology, for the banner.
inline const char* topology_name(input::TopologyKind k)
{
    switch (k){
    case input::TopologyKind::Global: return "global";
    case input::TopologyKind::Ring:   return "ring";
    case input::TopologyKind::Knn:    return "knn";
    case input::TopologyKind::Lsh:    return "lsh";
    }
    return "unknown";
}

/// @brief Reads the input file, validates it, runs the search and reports.
inline int run(const char* label, const char* input_path)
{
    input::Parameters p;

    try{
        p = input::read_parameters(input_path);
    }
    catch (const std::exception& error){
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }

    for (const auto& problem : input::validate(p)){
        std::cerr << "error: " << problem << '\n';
        return EXIT_FAILURE;
    }

    set_long_real_precision(p.convergence.mpfr_prec);

    std::cout << std::setprecision(16);
    std::cout << label << '\n'
              << "  input       " << input_path << '\n'
              << "  strategy    " << strategy_name(p.strategy) << '\n'
              << "  topology    " << topology_name(p.topology.kind) << '\n'
              << "  integrator  " << p.objective.poincare.integrator << '\n'
              << "  swap_tol    " << p.convergence.swap_tol
              << "   (" << p.convergence.mpfr_prec << "-bit exact phase)\n"
              << "  threads     " << p.execution.threads << '\n';

    const auto started = std::chrono::steady_clock::now();

    driver::Minima found;

    try{
        found = driver::run(p, std::cout);
    }
    catch (const std::exception& error){
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }

    const auto elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - started).count();

    const size_t stride = objective_stride(p);
    const size_t count  = found.size() / stride;

    std::cout << "\nwall clock  " << elapsed << " s\n"
              << count << (count == 1 ? " minimum" : " minima") << " found\n";

    for (size_t m = 0; m < count; m++){
        std::cout << "  orbit " << (m + 1) << '\n';
        for (size_t s = 0; s < stride; s++){
            const auto& point = found[m*stride + s];
            std::cout << "    x = " << point[0] << "   px = " << point[1] << '\n';
        }
    }

    return found.empty() ? EXIT_FAILURE : EXIT_SUCCESS;
}

} // namespace demo

#endif // PSO_DEMOS_DEMO_MAIN_HPP
