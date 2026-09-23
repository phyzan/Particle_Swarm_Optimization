// Analogue of local_classic_demo.cpp -- LSH neighbourhoods.
//
// Parameters live in demos/local_demo.toml, never in this file.

#include <demos/demo_main.hpp>

#ifndef PSO_DEMO_INPUT
#define PSO_DEMO_INPUT "demos/local_demo.toml"
#endif

int main()
{
    return demo::run("local_demo", PSO_DEMO_INPUT);
}
