#ifndef __DEFI__
#define __DEFI__

#include "../enhanced.hpp"

template <typename Type, typename Type_Arr, typename Type_Vec, typename Type_Empty>
class Deflection_Internal : public Enhanced<Type, Type_Arr, Type_Vec, Type_Empty>
{
  public:
    Deflection_Internal(const pso_params<Type, Type_Arr> &p, const enhanced_params<Type> &ep, std::string var_type,
                        Type l, Type constant, int precision, bool constriction = false, bool LLSH = false,
                        bool repulsion = false, std::ostream *output = &(std::cout))
        : Enhanced<Type, Type_Arr, Type_Vec, Type_Empty>(p, var_type, precision, constriction, LLSH, repulsion, output)
    {
        this->l = l;
        this->ep = ep;
        this->constant = constant;
    }

    bool parameter_check() override
    {
        // A function that checks if the parameters provided are
        // logically correct.

        if (!Enhanced<Type, Type_Arr, Type_Vec, Type_Empty>::parameter_check())
        {
            return false;
        }

        if (this->l <= 0)
        {
            std::cout << "~> Error: The deflection parameter l should be greater than 0." << std::endl;

            return false;
        }

        return true;
    }

    void print_params() override
    {
        // A functions that ouputs all the parameters used in the execution.

        Enhanced<Type, Type_Arr, Type_Vec, Type_Empty>::print_params();

        (*this->output) << "/- Deflection Parameters -" << std::endl;
        (*this->output) << "|- l        : " << this->l << std::endl;
        (*this->output) << "|- Constant : " << this->constant << std::endl;
        (*this->output) << "\\-------------------------" << std::endl;
    }

    Type_Arr fit(double swap_point, Type_Arr result)
    {
        bool success = false;

        this->result = result;

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

            return this->result;
        }

        this->fpopul = Type_Vec(this->p.popsize);

        this->initialize_arrays();
        this->evaluate_swarm("initial");

        if (!this->check_initial_particles())
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

        // Reset success to false in order to keep the program running,
        // since deflection can find multiple orbits.
        success = false;

        if ((swap_point - double(this->p.gm) <= double(this->p.gm) * 10))
        {
            this->obj_print_result(this->result);
        }

        return this->result;
    }

  private:
    Type l;
    Type constant; // The value of the minimum that the algorithm converges to.

    Type_Vec Pi(const Type_Arr& popul)
    {
        // Please refer to the publication for this function.

        Type_Vec Ones_vec = Type_Vec::Ones(this->p.popsize);
        Type_Vec Pi_result = Type_Vec::Ones(this->p.popsize);

        if (this->result.size() == 0)
        {
            return Pi_result;
        }
        else
        {
            for (int i = 0; i < this->result.cols(); i++)
            {
                Type_Arr curr_minimum(this->p.dim, this->p.popsize);

                curr_minimum = this->result.col(i).replicate(1, this->p.popsize);
                Pi_result = Pi_result.array() * (Ones_vec.array() / this->T(popul, curr_minimum).array()).array();
            }

            return Pi_result;
        }
    }

    Type_Vec T(const Type_Arr& x, const Type_Arr& x_star)
    {
        // Please refer to the publication for this function.

        return ((this->l * (x - x_star).colwise().norm()).array()).tanh();
    }

    Type_Vec obj_function(const Type_Arr& particles) override
    {
        return ((this->Obj_F->call(particles).array() + this->constant) * this->Pi(particles).array());
    }
};

#endif