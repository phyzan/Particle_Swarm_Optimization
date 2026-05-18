#ifndef __GRID__
#define __GRID__

#include "../pso/enhanced/local_classic/local_classic.hpp"
#include "../pso/global_classic/global_classic.hpp"

struct grid_params
{
    int resolution = 1;
    int batch_size = 1;
    int starting_subspace = 0;
    std::string dest = "test";
    std::string type = "Global_Classic";

    obj_params<MP_REAL> objp;
    enhanced_params<MP_REAL> ep;
    pso_params<MP_REAL, MPFR_ARR> p;
};

void init_folder(std::string &name);

MPFR_ARR *create_grid_bounds(const MPFR_ARR& original_bounds, int res, int &grid_size);

// Note: gp is passed by value intentionally as it's modified internally
MPFR_ARR grid_search(grid_params gp, double swap_point, int precision = 64, const std::string& Obj_F = "Poincare",
                     bool constriction = false);

#endif