#include "../src/pso/enhanced/deflection/deflection.hpp"

int main()
{
    MPFR_ARR bounds(2, 2);

    obj_params<MP_REAL> objp;
    enhanced_params<MP_REAL> ep;
    pso_params<MP_REAL, MPFR_ARR> p;

    bounds(0, 0) = -2.5;
    bounds(0, 1) = 2.5;
    bounds(1, 0) = -2.5;
    bounds(1, 1) = 2.5;

    p.gm = 0;
    p.dim = 2;
    p.min_w = 0.01;
    p.max_w = 0.5;
    p.c1 = 2;
    p.c2 = 1.7;
    p.max_it = 15000;
    p.bounds = bounds;
    p.popsize = 20;
    p.err_goal = 1e-8;

    ep.rep_radius = 0.25;
    ep.rep_rho = 3;

    objp.pc.p = 1;
    objp.pc.ene = 17;
    objp.pc.xpoin = -1.8019693;
    objp.pc.threads = 15;
    objp.pc.err_goal = 1e-8;

    Deflection *pso = Deflection_Init(p, ep, 1e-5, 50, 0, 100, false, false, true);

    if (pso == nullptr)
    {
        return EXIT_SUCCESS;
    }

    pso->print_params();
    pso->set_Obj_F("Poincare", objp);
    pso->fit(3);

    delete pso;
    pso = nullptr;

    mpfr_free_cache();

    return EXIT_SUCCESS;
}
