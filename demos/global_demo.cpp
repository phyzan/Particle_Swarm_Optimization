// Analogue of global_classic_demo.cpp -- one global attractor.
//
// Parameters live in demos/global_demo.toml, never in this file.

#include <demos/demo_main.hpp>

#ifndef PSO_DEMO_INPUT
#define PSO_DEMO_INPUT "demos/global_demo.toml"
#endif

int main()
{
    return demo::run("global_demo", PSO_DEMO_INPUT);
}
