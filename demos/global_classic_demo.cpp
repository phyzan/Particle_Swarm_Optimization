#include "../src/pso/global_classic/global_classic.hpp"

int main()
{
    MPFR_ARR bounds(2, 2);

    obj_params<MP_REAL> objp;
    pso_params<MP_REAL, MPFR_ARR> p;

    bounds(0, 0) = -2.5;
    bounds(0, 1) = 2.5;
    bounds(1, 0) = -0.5;
    bounds(1, 1) = 0.5;

    p.gm = 0;
    p.dim = 2;
    p.max_w = 0.5;
    p.c1 = 2;
    p.c2 = 1.7;
    p.max_it = 15000;
    p.bounds = bounds;
    p.popsize = 20;
    p.err_goal = 1e-8;

    objp.pc.p = 1;
    objp.pc.ene = 17;
    objp.pc.xpoin = -1.8019693;
    objp.pc.threads = 15;
    objp.pc.err_goal = 1e-8;

    Global_Classic *pso = Global_Classic_Init(p, 1e-5);

    if (pso == nullptr)
    {
        return EXIT_SUCCESS;
    }

    pso->print_params();
    pso->set_Obj_F("Poincare", objp);
    pso->fit();

    delete pso;
    pso = nullptr;

    mpfr_free_cache();

    return EXIT_SUCCESS;
}
