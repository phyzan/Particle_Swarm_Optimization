#include "global_classic.hpp"

Global_Classic::Global_Classic(pso_params<MP_REAL, MPFR_ARR> &p, double swap_point, int precision, std::ostream *output)
{
    this->swap_point = swap_point;

    // Initialize the MPFR model.

    this->gclc_in_mpreal =
        new Global_Classic_Internal<MP_REAL, MPFR_ARR, MPFR_VEC, MPFR_EMPTY>(p, "mp_real", precision, output);

    if (this->parameter_check() == false)
    { // If parameter check fails, return.
        this->initialized = false;
        this->gclc_in_double = nullptr;

        return;
    }
    else
    {
        this->initialized = true;
    }

    // Initialize the Double model.

    this->gclc_in_double = new Global_Classic_Internal<double, DOUBLE_ARR, DOUBLE_VEC, DOUBLE_EMPTY>(
        pso_params_to_double(p), "double", precision, output);
}

Global_Classic::~Global_Classic()
{
    if (this->gclc_in_double)
    {
        delete this->gclc_in_double;
        this->gclc_in_double = nullptr;
    }

    if (this->gclc_in_mpreal)
    {
        delete this->gclc_in_mpreal;
        this->gclc_in_mpreal = nullptr;
    }
}

void Global_Classic::set_Obj_F(const std::string &objective_func, const obj_params<MP_REAL> &params)
{
    // Set the objective functions for both models.

    this->gclc_in_mpreal->set_Obj_F(objective_func, params);
    this->gclc_in_double->set_Obj_F(objective_func, obj_params_to_double(params));
}

MPFR_ARR Global_Classic::fit()
{
    DOUBLE_ARR double_result;

    // Run the algorithm up until the swap point.
    double_result = this->gclc_in_double->fit(this->swap_point);

    if (double_result.size() > 0)
    {
        // If a minimum is found, return it and exit.
        MPFR_ARR mp_result;
        return double_to_mpfr_ARR(mp_result, double_result);
    }

    // Import the variables from the double model to the mpreal model.
    static thread_local variables v;
    this->gclc_in_mpreal->var_import(this->gclc_in_double->var_export(v));

    // Run from the swap point until the algorithm converges or
    // reaches maximum iterations, return the minimum (if one is found) and exit.
    return this->gclc_in_mpreal->fit(this->swap_point);
}

Global_Classic *Global_Classic_Init(pso_params<MP_REAL, MPFR_ARR> &p, double swap_point, int precision,
                                    std::ostream *output)
{
    // A function that initializes the Global Classic algorithm and checks that the
    // parameters are logically correct.

    Global_Classic *pso = new Global_Classic(p, swap_point, precision, output);

    if (pso->is_initialized() == false)
    {
        delete pso;
        pso = nullptr;
    }

    return pso;
}
