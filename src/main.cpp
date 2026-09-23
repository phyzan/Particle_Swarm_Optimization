#include <src/driver/dispatch.hpp>
#include <src/input/initialize.hpp>
#include <src/input/validate.hpp>
#include <src/objective/objective.hpp>

#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {

enum class Exit : int{
    Converged  = 0,
    NotFound   = 1,
    BadInput   = 2
};

/// @brief Prints every minimum found, grouped into orbits.
void report(const input::Parameters& p, const driver::Minima& found, std::ostream& out)
{
    const size_t stride = objective_stride(p);
    const size_t count  = found.size() / stride;

    out << "\n" << count << (count == 1 ? " minimum" : " minima") << " found\n";

    for (size_t m = 0; m < count; m++){
        out << "  orbit " << (m + 1) << '\n';
        for (size_t s = 0; s < stride; s++){
            const auto& point = found[m*stride + s];
            out << "    x = " << point[0] << "   px = " << point[1] << '\n';
        }
    }
}

} // namespace


/// @brief Reads the input file, validates it, runs the search and reports.
int main(int argc, char** argv)
{
    const std::string path = (argc > 1) ? argv[1] : "input.toml";

    input::Parameters params;

    try{
        params = input::read_parameters(path);
    }
    catch (const std::exception& error){
        std::cerr << "error: " << error.what() << '\n';
        return int(Exit::BadInput);
    }

    const auto problems = input::validate(params);

    if (!problems.empty()){
        std::cerr << "error: invalid parameters in " << path << ":\n";
        for (const auto& problem : problems){
            std::cerr << "  " << problem << '\n';
        }
        return int(Exit::BadInput);
    }

    set_long_real_precision(params.convergence.mpfr_prec);

    std::cout << std::setprecision(16);

    driver::Minima found;

    try{
        found = driver::run(params, std::cout);
    }
    catch (const std::exception& error){
        std::cerr << "error: " << error.what() << '\n';
        return int(Exit::NotFound);
    }

    report(params, found, std::cout);

    return found.empty() ? int(Exit::NotFound) : int(Exit::Converged);
}
