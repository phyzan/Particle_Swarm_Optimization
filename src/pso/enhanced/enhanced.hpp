#ifndef __ENH__
#define __ENH__

#include "../../lsh/lsh.hpp"
#include "../pso.hpp"

template <typename Type> struct enhanced_params
{
    // A structure that contains variables used by the enhanced algorithms
    // of PSO.
    // Please refer to the publication for further information on the
    // functionality of each variable.

    int lsh_k = 5;
    int lsh_w = 3;
    int lsh_L = 30;
    Type con_k = 1;
    Type rep_rho = 0.2;
    Type rep_radius = 0.1;
    Type lsh_radius = 0.2;
};

template <typename Type> enhanced_params<MP_REAL> enhanced_params_to_mpfr(enhanced_params<Type> params)
{
    enhanced_params<MP_REAL> mpfr_params;

    mpfr_params.lsh_k = params.lsh_k;
    mpfr_params.lsh_w = params.lsh_w;
    mpfr_params.lsh_L = params.lsh_L;
    mpfr_params.con_k = MP_REAL(params.con_k);
    mpfr_params.rep_rho = MP_REAL(params.rep_rho);
    mpfr_params.rep_radius = MP_REAL(params.rep_radius);
    mpfr_params.lsh_radius = MP_REAL(params.lsh_radius);

    return mpfr_params;
}

template <typename Type> enhanced_params<double> enhanced_params_to_double(enhanced_params<Type> params)
{
    enhanced_params<double> double_params;

    double_params.lsh_k = params.lsh_k;
    double_params.lsh_w = params.lsh_w;
    double_params.lsh_L = params.lsh_L;
    double_params.con_k = double(params.con_k);
    double_params.rep_rho = double(params.rep_rho);
    double_params.rep_radius = double(params.rep_radius);
    double_params.lsh_radius = double(params.lsh_radius);

    return double_params;
}

template <typename Type, typename Type_Arr, typename Type_Vec, typename Type_Empty>
class Enhanced : public PSO<Type, Type_Arr, Type_Vec, Type_Empty>
{
  public:
    Enhanced(const pso_params<Type, Type_Arr> &p, std::string var_type, int precision, bool constriction, bool LLSH,
             bool repulsion, std::ostream *output)
        : PSO<Type, Type_Arr, Type_Vec, Type_Empty>(p, var_type, precision, output)
    {
        this->lsh_flag = LLSH;
        this->rep_flag = repulsion;
        this->con_flag = constriction;

        if (this->con_flag == true)
        {
            this->constriction_init();
        }
    };

    bool parameter_check() override
    {
        // A function that checks if the parameters provided are
        // logically correct.

        if (PSO<Type, Type_Arr, Type_Vec, Type_Empty>::parameter_check() == false)
        {
            return false;
        }

        if (this->con_flag == true)
        {
            if (this->ep.con_k == 0)
            {
                std::cout << "~> Error: The constriction coefficient (k) should not be equal to 0." << std::endl;

                return false;
            }
            else if (this->constriction_init(true) == false)
            {

                return false;
            }
        }

        if (this->lsh_flag == true)
        {
            if (this->ep.lsh_k < 4)
            {
                std::cout << "~> Error: In order for LSH to work reliably, the number of  " << std::endl;
                std::cout << "function h() in each hash function g() should be at least 4." << std::endl;

                return false;
            }
            else if (this->ep.lsh_w < 2)
            {
                std::cout << "~> Error: W coefficient in LSH should be greater or equal to 2." << std::endl;

                return false;
            }
            else if (this->ep.lsh_L < 5)
            {
                std::cout << "~> Error: For LSH to work reliably there should exist at least 5 hash tables (L)."
                          << std::endl;

                return false;
            }
            else if (this->ep.lsh_radius <= 0 || this->ep.lsh_radius > 1)
            {
                std::cout << "~> Error: The LSH effective radius should be greater than 0 and smaller or equal to 1."
                          << std::endl;

                return false;
            }
        }

        if (this->rep_flag == true)
        {
            if (this->ep.rep_rho == 0)
            {
                std::cout << "~> Error: Rho should not be equal to 0." << std::endl;

                return false;
            }
            else if (this->ep.rep_radius <= 0)
            {
                std::cout << "~> Error: The repulsion radius should be greater than 0." << std::endl;

                return false;
            }
        }

        return true;
    }

    void print_params() override
    {
        // A functions that ouputs the PSO Classic and Enhanced parameters used in the execution.

        PSO<Type, Type_Arr, Type_Vec, Type_Empty>::print_params();

        (*this->output) << "/-- Enhanced Parameters --" << std::endl;
        (*this->output) << "|- lsh k            : " << this->ep.lsh_k << std::endl;
        (*this->output) << "|- lsh L            : " << this->ep.lsh_L << std::endl;
        (*this->output) << "|- lsh radius       : " << this->ep.lsh_radius << std::endl;
        (*this->output) << "|- constriction k   : " << this->ep.con_k << std::endl;
        (*this->output) << "|- repulsion rho    : " << this->ep.rep_rho << std::endl;
        (*this->output) << "|- repulsion radius : " << this->ep.rep_radius << std::endl;
        (*this->output) << "\\-------------------------" << std::endl;
    }

  protected:
    bool con_flag;  // Constriction flag.
    bool lsh_flag;  // LSH flag.
    bool rep_flag;  // Repulsion flag.
    Type con_coeff; // Constriction coefficient.

    enhanced_params<Type> ep;

    bool constriction_init(bool is_check = false)
    {
        // A function that initializes the constriction coefficient.
        // If conditions are not met, return error.

        if (this->p.c1 + this->p.c2 <= 4)
        {
            if (is_check == true)
            {
                (*this->output)
                    << "~> Error: It must be true that c1 + c2 > 4 for the constriction coefficient to be calculated."
                    << std::endl;
            }
            return false;
        }
        else
        {
            Type fi = this->p.c1 + this->p.c2;

            this->con_coeff = lmath::abs((2 * this->ep.con_k) / (2 - fi - lmath::sqrt(lmath::pow(fi, 2) - 4 * fi)));

            return true;
        }
    }

    void repulsion()
    {
        // The repulsion function that can be used alongside deflection.
        // Please refer to the publication for more information.

        if (this->result.size() > 0)
        {
            for (int i = 0; i < int(this->result.cols()); i++)
            {
                Type_Vec d(this->p.popsize);
                Type_Arr minimum_array(this->p.dim, this->p.popsize);

                minimum_array = this->result.col(i).replicate(1, this->p.popsize);
                d = (minimum_array - this->popul).colwise().norm();

                for (int j = 0; j < int(this->p.popsize); j++)
                {
                    if (d(j) <= this->ep.rep_radius)
                    { // If a particle is in the repulsion radius,
                        Type_Vec z(this->p.dim);

                        z = ((this->popul.col(j) - minimum_array.col(j)) /
                             d(j)); // find the unitary vector z
                                    // with direction opposite of the found minimum,
                        this->popul.col(j) = this->popul.col(j) +
                                             this->ep.rep_rho * z.array(); // multiply it with the repulsion strength
                                                                           // (rho) and add it to the particle.
                        if (this->bestpos.size() > 0)
                        {
                            this->bestpos.col(j) = this->popul.col(j);
                        }
                        if (this->fbestpos.size() > 0)
                        {
                            this->fbestpos(j) = lmath::get_infinity<Type>();
                        }
                    }
                }
            }
        }
    }

    void constriction()
    {
        // Apply velocity constriction.

        this->vel = this->con_coeff * this->vel.array();
    }

    void swarm_update() override
    {
        this->popul = this->popul + this->vel;
        this->repulsion(); // Repulsion is applied after the particle update.
    }

    void velocity_update(Type w) override
    {
        int dim = this->p.dim;
        int popsize = this->p.popsize;

        Type_Arr A(dim, popsize);
        Type_Arr N(dim, popsize);
        Type_Arr R1(dim, popsize);
        Type_Arr R2(dim, popsize);

        if (this->lsh_flag == true)
        {
            // If the LSH flag is true, use LSH to calculate the new particle velocity.

            LSH<Type, Type_Arr, Type_Vec> *lsh =
                new LSH<Type, Type_Arr, Type_Vec>(this->ep.lsh_k, this->ep.lsh_w, this->ep.lsh_L, this->p.err_goal,
                                                  this->ep.lsh_radius, this->p.bounds, this->popul);

            for (int i = 0; i < popsize; i++)
            {
                // Initialize the best distance of a particle to infinity,
                Type best_dist = lmath::get_infinity<Type>();

                // find its closest neighbours (their indeces) and
                Eigen::Vector<int, Eigen::Dynamic> particle_neighbours = lsh->find_idx(this->popul.col(i));

                // save the closest particle (that is not the particle itself)
                // in the N array.
                for (int j = 0; j < particle_neighbours.size(); j++)
                {
                    if (i == particle_neighbours(j) && particle_neighbours.size() != 1)
                    {
                        continue;
                    }

                    if (this->fpopul(particle_neighbours(j)) < best_dist)
                    {
                        N.col(i) = this->popul.col(particle_neighbours(j));
                        best_dist = this->fpopul(particle_neighbours(j));
                    }
                }
            }

            delete lsh;
            lsh = nullptr;
        }
        else
        {
            // If LSH is not set to true, proceed with the classic PSO algorithm.
            A = this->bestpos.col(this->g).replicate(1, this->bestpos.cols());

            for (int i = 0; i < this->bestpos.cols(); i++)
            {
                A.col(i) = this->bestpos.col(this->g);
            }
        }

        R1 = Type_Arr::NullaryExpr(dim, popsize, [&]() { return Type(uniform_real(0, 1)); });
        R2 = Type_Arr::NullaryExpr(dim, popsize, [&]() { return Type(uniform_real(0, 1)); });

        if (this->lsh_flag == true)
        {
            this->vel =
                w * this->vel + this->p.c1 * R1 * (this->bestpos - this->popul) + this->p.c2 * R2 * (N - this->popul);
        }
        else
        {
            this->vel =
                w * this->vel + this->p.c1 * R1 * (this->bestpos - this->popul) + this->p.c2 * R2 * (A - this->popul);
        }

        this->vel_clamp();

        if (this->con_flag == true)
        {
            this->constriction();
        }
    }

    void evaluate_swarm(const std::string &type = "") override
    {
        // Evaluate position.

        // Initialize best position matrices.
        if (type == "initial")
        {
            this->repulsion();
            this->fpopul = this->obj_function(this->popul);

            this->bestpos = this->popul;
            this->fbestpos = this->fpopul;
        }
        else
        {
            this->fpopul = this->obj_function(this->popul);

            for (int i = 0; i < this->fpopul.rows(); i++)
            {
                if (this->fpopul(i) < this->fbestpos(i))
                {
                    this->fbestpos(i) = this->fpopul(i);
                    this->bestpos.col(i) = this->popul.col(i);
                }
            }
        }

        // Find best particle in the population.
        this->fbestpos.minCoeff(&(this->g));
    }
};

#endif