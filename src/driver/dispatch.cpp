#include <src/driver/dispatch.hpp>

#include <src/driver/grid.hpp>
#include <src/driver/search.hpp>

namespace driver{

/// @brief Runs the strategy the parameters select.
Minima run(const input::Parameters& params, std::ostream& out)
{
    switch (params.strategy){

    case input::SearchStrategy::Deflection:
        return run_deflection(params, out);

    case input::SearchStrategy::Grid:
        return run_grid(params, out);

    case input::SearchStrategy::Single:
        break;
    }

    return run_single(params, out);
}

} // namespace driver
