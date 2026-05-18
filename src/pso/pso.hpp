#ifndef __PSO__
#define __PSO__

#include <chrono>

#include "../functions/functions.hpp"
#include "../objective_functions/objective_functions.hpp"


template <typename Type, typename Type_Arr> struct pso_params
{
    // A structure that contains variables used by all algorithms
    // of PSO.
    // Please refer to the publication for further information on
    // the functionality of each variable.

    int dim = 2;
    int max_it = 2000;
    int popsize = 20;
    Type c1 = 0.5;
    Type c2 = 0.5;
    Type gm = 0;
    Type max_w = 0.5;
    Type min_w = 0.01;
    Type err_goal = 1e-5;
    double swap_point;
    Type_Arr bounds;
};

struct variables
{
    // A structure containing essential variables and arrays
    // needed to resume the execution in the second (mp_real)
    // model, after the swap point in reached in the double
    // model.

    int iter;
    MPFR_ARR vel;
    MPFR_ARR popul;
    MPFR_VEC fpopul;
    MPFR_ARR result;
    MPFR_ARR bestpos;
    MPFR_VEC fbestpos;
    Eigen::Index g;
};

template <typename Type, typename Type_Arr>
pso_params<MP_REAL, MPFR_ARR> pso_params_to_mpfr(const pso_params<Type, Type_Arr>& params)
{
    pso_params<MP_REAL, MPFR_ARR> mpfr_params;

    mpfr_params.dim = params.dim;
    mpfr_params.max_it = params.max_it;
    mpfr_params.popsize = params.popsize;
    mpfr_params.c1 = MP_REAL(params.c1);
    mpfr_params.c2 = MP_REAL(params.c2);
    mpfr_params.gm = MP_REAL(params.gm);
    mpfr_params.max_w = MP_REAL(params.max_w);
    mpfr_params.min_w = MP_REAL(params.min_w);
    mpfr_params.err_goal = MP_REAL(params.err_goal);
    mpfr_params.bounds = double_to_mpfr_ARR(params.bounds);

    return mpfr_params;
}

template <typename Type, typename Type_Arr>
pso_params<double, DOUBLE_ARR> pso_params_to_double(const pso_params<Type, Type_Arr>& params)
{
    pso_params<double, DOUBLE_ARR> double_params;

    double_params.dim = params.dim;
    double_params.max_it = params.max_it;
    double_params.popsize = params.popsize;
    double_params.c1 = double(params.c1);
    double_params.c2 = double(params.c2);
    double_params.gm = double(params.gm);
    double_params.max_w = double(params.max_w);
    double_params.min_w = double(params.min_w);
    double_params.err_goal = double(params.err_goal);
    mpfr_to_double_ARR(double_params.bounds, params.bounds);

    return double_params;
}

template <typename Type, typename Type_Arr, typename Type_Vec, typename Type_Empty> class PSO
{
  public:
    PSO(const pso_params<Type, Type_Arr> &p, const std::string& var_type, int precision, std::ostream *output)
    {
        this->g = -1;
        this->p = p;
        this->Obj_F = nullptr;
        this->result = Type_Arr(0, 0);
        this->output = output;
        this->imported = false;
        this->var_type = var_type;

        if (precision < 64)
        {
            (*this->output) << "Default MP-Real precision cannot be lower than 64." << std::endl;
            (*this->output) << "Default precision set to 64." << std::endl;

            precision = 64;
        }

        if (this->p.dim <= 0)
        {
            return;
        }

        // Sets default MPFR precision globally.
        lazy::set_default_mpreal_prec(precision);
        // mpfr::mpreal::set_default_prec(precision);

        // The velocity boundaries are arbitrarily set boundaries for the
        // velocity of the particles with the aim to prevent numerical explosions.
        // For each dimension, the maximum velocity of a particle can be half of the
        // distance of the min and max boundary of said dimension.
        // For example if for dimension 0, the minimum value is 5 and the maximum is 10,
        // the maximum velocity change that can occur in a single iteration for a particle
        // would be |(10-5)/2| = |2.5|, for dimension 0.
        this->vel_bounds = Type_Arr(this->p.dim, 2);

        for (int i = 0; i < this->p.dim; i++)
        {
            Type abs_vebound = abs((this->p.bounds(i, 1) - this->p.bounds(i, 0)) / 2);

            this->vel_bounds(i, 0) = -abs_vebound;
            this->vel_bounds(i, 1) = abs_vebound;
        }

        // Sets a global precision for the printed values.
        (*this->output) << std::setprecision(16);
    }

    virtual ~PSO()
    {
        if (this->Obj_F != nullptr)
        {
            delete this->Obj_F;
            this->Obj_F = nullptr;
        }
    }

    void var_import(const variables& v)
    {
        // A function that uses a variables structure to import
        // essential variables to continue the execution with MPFR variables.

        this->iter = v.iter;
        this->vel = v.vel;
        this->popul = v.popul;
        this->fpopul = v.fpopul;
        this->save_minima(v.result);
        this->bestpos = v.bestpos;
        this->fbestpos = v.fbestpos;
        this->g = v.g;

        this->imported = true;
    }

    void set_Obj_F(const std::string& objective_func, const obj_params<Type> &params)
    {
        this->Obj_F = new Objective_Functions<Type, Type_Arr, Type_Vec, Type_Empty>(this->output);

        this->Obj_F->init(params, objective_func);
    }

    variables var_export(variables& v)
    {
        // A function that creates a variables structure and exports
        // the essential variables.

        v.iter = this->iter;
        double_to_mpfr_ARR(v.vel, this->vel);
        double_to_mpfr_ARR(v.popul, this->popul);
        double_to_mpfr_VEC(v.fpopul, this->fpopul);
        double_to_mpfr_ARR(v.result, this->result);
        double_to_mpfr_ARR(v.bestpos, this->bestpos);
        double_to_mpfr_VEC(v.fbestpos, this->fbestpos);
        v.g = this->g;

        return v;
    }

    virtual bool parameter_check()
    {
        // A function that checks if the parameters provided are
        // logically correct.

        if (this->p.dim <= 0)
        {
            std::cout << "~> Error: Dimension cannot be 0 or a negative number." << std::endl;

            return false;
        }
        else if (this->p.max_it <= 0)
        {
            std::cout << "~> Error: The maximum iterations should be at least 1." << std::endl;

            return false;
        }
        else if (this->p.popsize <= 0)
        {
            std::cout << "~> Error: The population size cannot be 0 or a negative number." << std::endl;

            return false;
        }
        else if (this->p.err_goal <= 0)
        {
            std::cout << "~> Error: The error goal (precision) must be positive (err_goal > 0)." << std::endl;

            return false;
        }
        else if (this->p.max_w < 0 || this->p.min_w < 0)
        {
            std::cout << "~> Error: Both the min and max weights must be greater or equal to 0." << std::endl;

            return false;
        }
        else if (this->p.max_w <= this->p.min_w)
        {
            std::cout << "~> Error: The maximum weight value must be greater than the minimum (max_w > min_w)."
                      << std::endl;

            return false;
        }

        for (int i = 0; i < this->p.dim; i++)
        {
            if (this->p.bounds(i, 0) > this->p.bounds(i, 1))
            {
                std::cout << "~> Error: The upper bound must be greater or equal to the lower bound." << std::endl;

                return false;
            }
        }

        return true;
    };

    virtual void print_params()
    {
        // A function that ouputs the PSO Classic parameters used in the execution.

        (*this->output) << "/---- PSO  Parameters ----" << std::endl;
        (*this->output) << "|- Dimensions      : " << this->p.dim << std::endl;
        (*this->output) << "|- Max Iterations  : " << this->p.max_it << std::endl;
        (*this->output) << "|- Population size : " << this->p.popsize << std::endl;
        (*this->output) << "|- C1              : " << this->p.c1 << std::endl;
        (*this->output) << "|- C2              : " << this->p.c2 << std::endl;
        (*this->output) << "|- Global Minimum  : " << this->p.gm << std::endl;
        (*this->output) << "|- Max weight      : " << this->p.max_w << std::endl;
        (*this->output) << "|- Min weight      : " << this->p.min_w << std::endl;
        (*this->output) << "|- Error Goal      : " << this->p.err_goal << std::endl;
        (*this->output) << "\\-------------------------" << std::endl;

        (*this->output) << "\\-------- Bounds ---------" << std::endl;
        (*this->output) << this->p.bounds << std::endl;
        (*this->output) << "\\-------------------------" << std::endl;
    }

  protected:
    int iter;
    bool imported;
    MP_REAL avg_time;
    Type_Arr vel;
    Type_Arr popul;
    Type_Vec fpopul;
    Type_Arr result;
    Type_Arr bestpos;
    Type_Vec fbestpos;
    Type_Arr vel_bounds;
    std::string var_type;
    Eigen::Index g;
    std::ostream *output;

    pso_params<Type, Type_Arr> p;

    Objective_Functions<Type, Type_Arr, Type_Vec, Type_Empty> *Obj_F;

    bool check_initial_particles()
    {
        // Returns false if all particles' energies are out-of-bound.

        int inf_count = 0;

        for (int i = 0; i < this->fpopul.size(); i++)
        {
            if (this->fpopul(i) == std::numeric_limits<Type>::infinity())
            {
                inf_count++;
            }
        }

        return (inf_count != this->fpopul.size());
    }

    void vel_clamp()
    {
        // A function that prevents the explosion of particle values
        // using the velocity boundaries as explained in the constructor.
        // Clamp in-place using vectorized operations - more cache-friendly.

        for (int j = 0; j < this->p.popsize; j++)
        {
            // Column-wise access is cache-friendly for column-major Eigen arrays
            for (int i = 0; i < this->p.dim; i++)
            {
                if (this->vel(i, j) < this->vel_bounds(i, 0)) {
                    this->vel(i, j) = this->vel_bounds(i, 0);
                } else if (this->vel(i, j) > this->vel_bounds(i, 1)) {
                    this->vel(i, j) = this->vel_bounds(i, 1);
                }
            }
        }
    }

    void initialize_arrays()
    {
        if (this->imported)
        {
            return;
        }

        Type_Empty empty;

        // Initialize velocity and population randomly with the U(0,1) distribution.
        this->vel = Type_Arr::NullaryExpr(this->p.dim, this->p.popsize, [&]() { return Type(uniform_real(0, 1)); });
        this->popul = Type_Arr::NullaryExpr(this->p.dim, this->p.popsize, [&]() { return Type(uniform_real(0, 1)); });

        // Apply dimensional limits for the population.
        // Column-wise access is cache-friendly for column-major Eigen arrays
        for (int j = 0; j < this->p.popsize; j++)
        {
            for (int i = 0; i < this->p.dim; i++)
            {
                this->popul(i, j) = this->popul(i, j) * (this->p.bounds(i, 1) - this->p.bounds(i, 0)) + this->p.bounds(i, 0);
            }
        }

        this->fpopul = empty;
        this->bestpos = empty;
        this->fbestpos = empty;
    }

    void save_minima(const Type_Arr& min)
    {
        // Expands the result array to save a new minimum that is found.
        // Uses conservativeResize to avoid creating temporary array.

        if (this->result.size() == 0)
        {
            this->result = min;
        }
        else
        {
            Eigen::Index old_cols = this->result.cols();
            // conservativeResize preserves existing data
            this->result.conservativeResize(this->p.dim, old_cols + min.cols());
            this->result.rightCols(min.cols()) = min;
        }
    }

    void check_stop_criterion(bool &success)
    {
        if (abs(this->fbestpos(g) - this->p.gm) <= this->p.err_goal)
        {
            success = true;
        }
    }

    void swarm_evolution(bool &success, double swap_point)
    {
        if (!(this->imported))
        {
            this->iter = 0;
        }

        Type w = this->p.max_w;
        Type weveryit = floor(0.75 * this->p.max_it);
        Type inertdec = (this->p.max_w - this->p.min_w) / weveryit;

        // Swarm evolution loop.
        while ((!success) && (this->iter < this->p.max_it))
        {
            auto start = std::chrono::high_resolution_clock::now();

            if (swap_point > this->fbestpos(g))
            {
                return;
            }

            this->iter++;

            // Update the value of the inertia weight w.
            if (this->iter <= weveryit)
            {
                w = this->p.max_w - (this->iter - 1) * inertdec;
            }

            this->velocity_update(w);
            this->swarm_update();
            this->evaluate_swarm();
            this->check_stop_criterion(success);

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

            this->avg_time += duration.count() / (MP_REAL)1000000000.0;

            if (success || this->iter >= this->p.max_it)
            {
                this->avg_time /= this->iter;
            }

            // (*this->output) << "|----------------- " << std::endl;
            // (*this->output) << "|- Iteration     : " << this->iter << std::endl;
            // (*this->output) << "|- Error         : " << this->fbestpos(g) << std::endl;
            // (*this->output) << "|- Duration(s)   : " << duration.count() / 1000000000.0 << std::endl;
            // (*this->output) << "|- Best Particle : " << this->bestpos.col(g).transpose() << std::endl;
            // (*this->output) << "|----------------- " << std::endl;
        }

        this->imported = false;
    }

    virtual void swarm_update()
    {
        // The default swarm update function.

        this->popul = this->popul + this->vel;
    }

    virtual void velocity_update(Type w)
    {
        // The default velocity update function.

        int dim = this->p.dim;
        int popsize = this->p.popsize;

        Type_Arr A(dim, popsize);
        Type_Arr R1(dim, popsize);
        Type_Arr R2(dim, popsize);

        // Use Eigen's replicate - no need for the redundant loop
        A = bestpos.col(g).replicate(1, popsize);

        R1 = Type_Arr::NullaryExpr(dim, popsize, [&]() { return Type(uniform_real(0, 1)); });
        R2 = Type_Arr::NullaryExpr(dim, popsize, [&]() { return Type(uniform_real(0, 1)); });
        vel = w * vel + this->p.c1 * R1 * (bestpos - popul) + this->p.c2 * R2 * (A - popul);

        this->vel_clamp();
    }

    virtual void evaluate_swarm(const std::string &type = "")
    {
        // Evaluate position.

        this->fpopul = this->obj_function(this->popul);

        // Initialize best position matrices.
        if (type == "initial" && !this->imported)
        {
            this->bestpos = this->popul;
            this->fbestpos = this->fpopul;
        }
        else
        {
            for (int i = 0; i < fpopul.rows(); i++)
            {
                if (this->fpopul(i) < this->fbestpos(i))
                {
                    this->fbestpos(i) = this->fpopul(i);
                    this->bestpos.col(i) = this->popul.col(i);
                }
            }
        }

        // Find the best particle in the population.
        this->fbestpos.minCoeff(&g);
    }

    inline void obj_print_result(const Type_Arr& res)
    {
        this->Obj_F->print_result(res);
    }

    inline Type_Arr obj_calculate(const Type_Vec& particle)
    {
        return this->Obj_F->calculate(particle);
    }

    inline virtual Type_Vec obj_function(const Type_Arr& population)
    {
        return this->Obj_F->call(population).array();
    }
};

#endif