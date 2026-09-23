// Analogue of deflection_demo.cpp -- repeated runs, poles at known minima.
//
// Parameters live in demos/deflection_demo.toml, never in this file.

#include <demos/demo_main.hpp>

#ifndef PSO_DEMO_INPUT
#define PSO_DEMO_INPUT "demos/deflection_demo.toml"
#endif

int main()
{
    return demo::run("deflection_demo", PSO_DEMO_INPUT);
}
