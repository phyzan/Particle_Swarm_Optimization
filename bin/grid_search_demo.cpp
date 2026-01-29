#include "../src/grid_search/grid_search.hpp"

int main()
{
    MPFR_ARR bounds(2, 2);

    grid_params gp;
    obj_params<MP_REAL> objp;
    enhanced_params<MP_REAL> ep;
    pso_params<MP_REAL, MPFR_ARR> p;

    bounds(0, 0) = -2.5;
    bounds(0, 1) = 2.5;
    bounds(1, 0) = -2.5;
    bounds(1, 1) = 2.5;

    p.gm = 0;
    p.max_w = 0.5;
    p.min_w = 0.01;
    p.c1 = 2;
    p.c2 = 1.7;
    p.max_it = 15000;
    p.bounds = bounds;
    p.popsize = 20;
    p.err_goal = 1e-8;

    objp.pc.p = 2;
    objp.pc.ene = 17;
    objp.pc.xpoin = -1.8019693;
    objp.pc.threads = 15;
    objp.pc.err_goal = 1e-8;

    gp.p = p;
    gp.objp = objp;
    gp.dest = "p2_1024";
    gp.type = "Global_Classic";
    gp.batch_size = 2;
    gp.resolution = 4;

    Poincare<MP_REAL, MPFR_ARR, MPFR_VEC> *dummy_poincare = new Poincare<MP_REAL, MPFR_ARR, MPFR_VEC>(objp.pc);

    dummy_poincare->print_results(grid_search(gp, 1e-5, 100, "Poincare"));

    delete dummy_poincare;
    dummy_poincare = nullptr;

    mpfr_free_cache();

    return EXIT_SUCCESS;
}
