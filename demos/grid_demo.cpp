// Analogue of grid_search_demo.cpp -- one swarm per subspace.
//
// Parameters live in demos/grid_demo.toml, never in this file.

#include <demos/demo_main.hpp>

#ifndef PSO_DEMO_INPUT
#define PSO_DEMO_INPUT "demos/grid_demo.toml"
#endif

int main()
{
    return demo::run("grid_demo", PSO_DEMO_INPUT);
}
