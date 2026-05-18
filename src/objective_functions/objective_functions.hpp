#ifndef __OBJF__
#define __OBJF__

#include "poincare_functions/poincare.hpp"
#include "poincare_functions/simple_poincare.hpp"


template <typename Type> struct obj_params
{
    poinc_params<Type> pc;
    // template_params<Type,Type_Arr,Type_Vec,Type_Empty> temp;
};

template <typename Type> obj_params<MP_REAL> obj_params_to_mpfr(const obj_params<Type>& params)
{
    obj_params<MP_REAL> mpfr_params;

    mpfr_params.pc.p = params.pc.p;
    mpfr_params.pc.threads = params.pc.threads;
    mpfr_params.pc.c1 = MP_REAL(params.pc.c1);
    mpfr_params.pc.c2 = MP_REAL(params.pc.c2);
    mpfr_params.pc.c3 = MP_REAL(params.pc.c3);
    mpfr_params.pc.dt = MP_REAL(params.pc.dt);
    mpfr_params.pc.ene = MP_REAL(params.pc.ene);
    mpfr_params.pc.xpoin = MP_REAL(params.pc.xpoin);
    mpfr_params.pc.err_goal = MP_REAL(params.pc.err_goal);

    /*
    Conversions for template parameters.
    */

    return mpfr_params;
}

template <typename Type> obj_params<double> obj_params_to_double(const obj_params<Type>& params)
{
    obj_params<double> double_params;

    double_params.pc.p = params.pc.p;
    double_params.pc.threads = params.pc.threads;
    double_params.pc.c1 = double(params.pc.c1);
    double_params.pc.c2 = double(params.pc.c2);
    double_params.pc.c3 = double(params.pc.c3);
    double_params.pc.dt = double(params.pc.dt);
    double_params.pc.ene = double(params.pc.ene);
    double_params.pc.xpoin = double(params.pc.xpoin);
    double_params.pc.err_goal = double(params.pc.err_goal);

    /*
    Conversions for template parameters.
    */

    return double_params;
}

template <typename Type, typename Type_Arr, typename Type_Vec, typename Type_Empty> class Objective_Functions
{
  public:
    Objective_Functions(std::ostream *output = &(std::cout))
    {
        this->output = output;

        this->poincare = nullptr;
        this->s_poincare = nullptr;
        /*
            this->template = nullptr;
        */
    }

    ~Objective_Functions()
    {
        if (this->poincare != nullptr)
        {
            delete this->poincare;
            this->poincare = nullptr;
        }

        if (this->s_poincare != nullptr)
        {
            delete this->s_poincare;
            this->s_poincare = nullptr;
        }
        /*
            if(this->template != nullptr){
                delete this->template;
                this->template = nullptr;
            }
        */
    }

    void print_result(const Type_Arr& result)
    {
        // A function to print the result in a custom way
        // tailored to the objective function that is in use.

        if (this->poincare != nullptr)
        {
            this->poincare->print_results(result);
        }

        if (this->s_poincare != nullptr)
        {
            this->s_poincare->print_results(result);
        }
        /*
            if(this->template != nullptr){
                this->template->print_results(result);
            }
        */
    }

    void init(const obj_params<Type> &params, const std::string &type)
    {
        if (type == "Poincare")
        {
            this->poincare = new Poincare<Type, Type_Arr, Type_Vec>(params.pc, this->output);
        }
        else if (type == "Simple_Poincare")
        {
            this->s_poincare = new Simple_Poincare<Type, Type_Arr, Type_Vec>(params.pc, this->output);
        }
        /*
        else if(type == "Template"){
            this->template = new Template<Type,Type_Arr,Type_Vec,Type_Empty>(params.temp);
        }
        */
    }

    Type_Vec call(const Type_Arr& popul)
    {
        // A function that calls the objective function in use.

        if (this->poincare != nullptr)
        {
            return this->poincare->objective_function(popul);
        }
        else if (this->s_poincare != nullptr)
        {
            return this->s_poincare->objective_function(popul);
        }
        /*
            else if(this->template != nullptr){
                this->template->objective_function(popul);
            }
        */
        else
        {
            Type_Empty empty;

            return empty;
        }
    }

    Type_Arr calculate(const Type_Vec& particle)
    {
        // A function that calculates the new position of only one particle of the population.

        if (this->poincare != nullptr)
        {
            Type_Vec q(4);

            return this->poincare->calculate(particle, q, false);
        }
        else if (this->s_poincare != nullptr)
        {
            Type_Vec q(4);

            return this->s_poincare->calculate(particle, q, false);
        }
        /*
            else if(this->template != nullptr){
                this->template->calculate(particle,false);
            }
        */
        else
        {
            Type_Empty empty;

            return empty;
        }
    }

  private:
    std::ostream *output;
    Poincare<Type, Type_Arr, Type_Vec> *poincare;
    Simple_Poincare<Type, Type_Arr, Type_Vec> *s_poincare;
    // Template                                *template;
};

#endif