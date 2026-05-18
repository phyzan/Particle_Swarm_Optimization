#ifndef __GCLI__
#define __GCLI__

#include "../pso.hpp"

template <typename Type, typename Type_Arr, typename Type_Vec, typename Type_Empty>
class Global_Classic_Internal : public PSO<Type, Type_Arr, Type_Vec, Type_Empty>
{
  public:
    Global_Classic_Internal(const pso_params<Type, Type_Arr> &p, std::string var_type, int precision = 64,
                            std::ostream *output = &(std::cout))
        : PSO<Type, Type_Arr, Type_Vec, Type_Empty>(p, var_type, precision, output) {};

    Type_Arr fit(double swap_point)
    {
        bool success = false;

        Type_Empty empty;
        this->result = empty;

        if (this->var_type == "mp_real")
        {
            // If the variable type is mp_real, set the swap point to the global minimum.
            // This allows the program to converge and exit.
            // If the variable type is double, then the program exits at the swap point
            // set by the user.

            swap_point = double(this->p.gm);
        }

        if (this->Obj_F == nullptr)
        {
            (*this->output) << "~> Error: Objective Function was not declared." << std::endl;

            return empty;
        }

        this->fpopul = Type_Vec(this->p.popsize);

        this->initialize_arrays();
        this->evaluate_swarm("initial");

        if (!(this->check_initial_particles()))
        {
            // If all the particles have out-of-bounds energy, return this error.

            (*this->output) << "Error: All initial particles are not suitable." << std::endl;

            return this->result;
        }

        this->swarm_evolution(success, swap_point);

        if (success)
        {
            // If an orbit is found, save it.

            this->save_minima(this->obj_calculate(this->bestpos.col(this->g)));

            (*this->output) << "|---------------------- " << std::endl;
            (*this->output) << "-     Orbit Found     - " << std::endl;
            (*this->output) << "-  Average Iteration  - " << std::endl;
            (*this->output) << "-    Time(seconds)    - " << std::endl;
            (*this->output) << "- " << this->avg_time << " - " << std::endl;
            (*this->output) << "|---------------------- " << std::endl;
        }

        if (swap_point - double(this->p.gm) <= double(this->p.gm) * 10)
        {
            this->obj_print_result(this->result);
        }

        return this->result;
    }
};

#endif