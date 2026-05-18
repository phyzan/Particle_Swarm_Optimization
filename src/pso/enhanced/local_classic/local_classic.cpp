#include "local_classic.hpp"

Local_Classic::Local_Classic(pso_params<MP_REAL, MPFR_ARR> &p, enhanced_params<MP_REAL> &ep, double swap_point,
                             int precision, bool constriction, std::ostream *output)
{
    this->swap_point = swap_point;

    // Initialize the MPFR model.

    this->lclc_in_mpreal = new Local_Classic_Internal<MP_REAL, MPFR_ARR, MPFR_VEC, MPFR_EMPTY>(
        p, ep, "mp_real", precision, constriction, output);

    if (this->parameter_check() == false)
    { // If parameter check fails, return.
        this->initialized = false;
        this->lclc_in_double = nullptr;

        return;
    }
    else
    {
        this->initialized = true;
    }

    pso_params<double, DOUBLE_ARR> dp; // A PSO parameter struct that contains all the variables, but in double.

    dp.dim = p.dim;
    dp.max_it = p.max_it;
    dp.popsize = p.popsize;
    dp.c1 = double(p.c1);
    dp.c2 = double(p.c2);
    dp.gm = double(p.gm);
    dp.max_w = double(p.max_w);
    dp.min_w = double(p.min_w);
    dp.err_goal = double(p.err_goal);

    DOUBLE_ARR dbounds(p.dim, 2);
    for (int i = 0; i < p.dim; i++)
    {
        dbounds(i, 0) = double(p.bounds(i, 0));
        dbounds(i, 1) = double(p.bounds(i, 1));
    }

    dp.bounds = dbounds;

    enhanced_params<double> dep; // An enhanced parameter struct that contains all the parameters, but in double.

    dep.lsh_k = ep.lsh_k;
    dep.lsh_w = ep.lsh_w;
    dep.lsh_L = ep.lsh_L;
    dep.con_k = double(ep.con_k);
    dep.rep_rho = double(ep.rep_rho);
    dep.rep_radius = double(ep.rep_radius);
    dep.lsh_radius = double(ep.lsh_radius);

    // Initialize the Double model.

    this->lclc_in_double = new Local_Classic_Internal<double, DOUBLE_ARR, DOUBLE_VEC, DOUBLE_EMPTY>(
        pso_params_to_double(p), enhanced_params_to_double(ep), "double", precision, constriction, output);
}

Local_Classic::~Local_Classic()
{
    if (this->lclc_in_double)
    {
        delete this->lclc_in_double;
        this->lclc_in_double = nullptr;
    }

    if (this->lclc_in_mpreal)
    {
        delete this->lclc_in_mpreal;
        this->lclc_in_mpreal = nullptr;
    }
}

void Local_Classic::set_Obj_F(const std::string &objective_func, const obj_params<MP_REAL> &params)
{
    // Set the objective functions for both models.

    this->lclc_in_mpreal->set_Obj_F(objective_func, params);
    this->lclc_in_double->set_Obj_F(objective_func, obj_params_to_double(params));
}

MPFR_ARR Local_Classic::fit()
{
    DOUBLE_ARR double_result;

    // Run the algorithm up until the swap point.
    double_result = this->lclc_in_double->fit(this->swap_point);

    if (double_result.size() > 0)
    {
        // If a minimum is found, return it and exit.
        MPFR_ARR mp_result;
        return double_to_mpfr_ARR(mp_result, double_result);
    }

    // Import the variables from the double model to the mpreal model.
    static thread_local variables v;
    this->lclc_in_mpreal->var_import(this->lclc_in_double->var_export(v));

    // Run from the swap point until the algorithm converges or
    // reaches maximum iterations, return the minimum (if one is found) and exit.
    return this->lclc_in_mpreal->fit(this->swap_point);
}

Local_Classic *Local_Classic_Init(pso_params<MP_REAL, MPFR_ARR> &p, enhanced_params<MP_REAL> &ep, double swap_point,
                                  int precision, bool constriction, std::ostream *output)
{
    // A function that initializes the Local Classic algorithm and checks that the
    // parameters are logically correct.

    Local_Classic *pso = new Local_Classic(p, ep, swap_point, precision, constriction, output);

    if (pso->is_initialized() == false)
    {
        delete pso;
        pso = nullptr;
    }

    return pso;
}
