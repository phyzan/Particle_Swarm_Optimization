#include "deflection.hpp"

Deflection::Deflection(pso_params<MP_REAL, MPFR_ARR> &p, enhanced_params<MP_REAL> &ep, double swap_point, MP_REAL &l,
                       MP_REAL &constant, int precision, bool constriction, bool LLSH, bool repulsion,
                       std::ostream *output)
{
    this->swap_point = swap_point;

    // Initialize the MPFR model.

    this->defl_in_mpreal = new Deflection_Internal<MP_REAL, MPFR_ARR, MPFR_VEC, MPFR_EMPTY>(
        p, ep, "mp_real", l, constant, precision, constriction, LLSH, repulsion, output);

    if (this->parameter_check() == false)
    { // If parameter check fails, return.
        this->initialized = false;
        this->defl_in_double = nullptr;

        return;
    }
    else
    {
        this->initialized = true;
    }

    // Initialize the Double model.

    this->defl_in_double = new Deflection_Internal<double, DOUBLE_ARR, DOUBLE_VEC, DOUBLE_EMPTY>(
        pso_params_to_double(p), enhanced_params_to_double(ep), "double", double(l), double(constant), precision,
        constriction, LLSH, repulsion, output);
}

Deflection::~Deflection()
{
    if (this->defl_in_double)
    {
        delete this->defl_in_double;
        this->defl_in_double = nullptr;
    }

    if (this->defl_in_mpreal)
    {
        delete this->defl_in_mpreal;
        this->defl_in_mpreal = nullptr;
    }
}

void Deflection::set_Obj_F(const std::string &objective_func, const obj_params<MP_REAL> &params)
{
    // Set the objective functions for both models.

    this->defl_in_mpreal->set_Obj_F(objective_func, params);
    this->defl_in_double->set_Obj_F(objective_func, obj_params_to_double(params));
}

MPFR_ARR Deflection::fit(int runs)
{
    MPFR_ARR mp_result;
    MPFR_EMPTY mp_empty;

    DOUBLE_ARR db_result;
    DOUBLE_EMPTY db_empty;

    mp_result = mp_empty;

    for (int i = 0; i < runs; i++)
    {
        // Initialize the results to empty arrays.
        if (mp_result.size() == 0)
        {
            db_result = db_empty;
        }
        else
        {
            db_result = mpfr_to_double_ARR(mp_result);
        }

        // Run the algorithm up until the swap point.
        db_result = this->defl_in_double->fit(this->swap_point, db_result);

        // If the maximum number (equal to the number of runs) of minima
        // is found, return them and exit.
        if (db_result.cols() == runs)
        {
            return double_to_mpfr_ARR(db_result);
        }

        // Import the variables from the double model to the mpreal model.
        this->defl_in_mpreal->var_import(this->defl_in_double->var_export());

        // Run from the swap point until the algorithm converges or
        // reaches maximum iterations.
        mp_result = double_to_mpfr_ARR(db_result);

        mp_result = this->defl_in_mpreal->fit(this->swap_point, mp_result);
    }

    if (db_result.cols() == runs)
    {
        return double_to_mpfr_ARR(db_result);
    }
    else
    {
        return mp_result;
    }
}

Deflection *Deflection_Init(pso_params<MP_REAL, MPFR_ARR> &p, enhanced_params<MP_REAL> &ep, double swap_point,
                            MP_REAL l, MP_REAL constant, int precision, bool constriction, bool LLSH, bool repulsion,
                            std::ostream *output)
{
    // A function that initializes the deflecion algorithm and checks that the
    // parameters are logically correct.

    Deflection *pso = new Deflection(p, ep, swap_point, l, constant, precision, constriction, LLSH, repulsion, output);

    if (pso->is_initialized() == false)
    {
        delete pso;
        pso = nullptr;
    }

    return pso;
}
