#ifndef __EXMP__
#define __EXMP__

/*
After creating the new class according to the following template,
complete the integration in the objective_functions.hpp file, as
instructed by the comments.
*/

#include "../../local_definitions.hpp"

template <typename Type, typename Type_Arr, typename Type_Vec, typename Type_Empty> struct template_params
{
    //--parameters--//
};

template <typename Type, typename Type_Arr, typename Type_Vec, typename Type_Empty> class Template
{
  public:
    Template(const template_params<Type, Type_Arr, Type_Vec, Type_Empty> &tp)
    {
        this->tp = tp;
    }

    ~Template()
    {
        // Free any customly allocated memory.
    }

    void print_results(Type_Arr results)
    {
        // Custom output based on the minimisation problem.
        // Could be left empty.
    }

    Type_Arr calculate(Type_Vec particle, bool first_only = true)
    {
        // Function that calculates the new position of one particle.
        // It can return multiple points in the space, or just the first one.
        // By default this should return the first point.
        // It should return the multiple points when it is called to calculate
        // the minima, for Template, when the minima points are saved.

        // Objective_function->calculate() -> calculate(_,false);
    }

    Type_Vec objective_function(Type_Arr population)
    {
        // This function should call this->calculate(_,true) for
        // every particle and measure the loss for each particle
        //(new_pos - population_pos).

        // This function should return an array with the calculated
        // loss for each function (1 * num_of_particles).
    }

  private:
    template_params<Type, Type_Arr, Type_Vec, Type_Empty> tp;
};

#endif