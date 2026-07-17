#ifndef __PPNC__
#define __PPNC__


#include <iostream>
#include <lazy/lazy.hpp>
#include <odepack/ode/Core/ObjectiveSolver.hpp>
#include <odepack/ode/Tools.hpp>
#include <odepack/odepack.hpp>
#include <lazy/mpfrLazy.hpp>
#include "../../local_definitions.hpp"

#ifdef __MAC__
#include <dispatch/dispatch.h>
#else
#include <mutex>
#include <omp.h>
#endif

template <typename Type> struct poinc_params
{
    // A structure containing variables that the user can set
    // for this Poincare implementation.

    int p = 1;
    int threads = 1;
    Type c1 = 5;
    Type c2 = 3;
    Type c3 = -0.3;
    Type dt = 0.001;
    Type ene = 17;
    Type xpoin = 0;
    Type err_goal = 1e-10;
};

template <typename Type, typename Type_Arr, typename Type_Vec> class Parent_Poincare;


template<typename T>
struct OdeSystem{

    T c1, c2, c3;

    NDSPAN_INLINE void Rhs(auto* out, const auto& t, const auto* q, const auto* args) const{

        const auto &x = q[0];
        const auto &y = q[1];
        const auto &px = q[2];
        const auto &py = q[3];

        auto vx = 2 * c1 * x - 4 * c3 * x * x * x + 12 * c3 * x * y * y;
        auto vy = 2 * c1 * y + c2 - 4 * c3 * y * y * y + 12 * c3 * y * x * x;

        out[0] = px;
        out[1] = py;
        out[2] = -vx;
        out[3] = -vy;
    }

};


template<typename T>
struct MyObjFunc{

    T point;

    inline T operator()(const T& t, const T* q, const T* /*args*/) const{
        return q[1] - point;
    }

};



template<typename T>
class MySolver : public ode::ObjectiveSolver<ode::RK45, T, 4, ode::SolverPolicy::Static, OdeSystem<T>, MySolver<T>, MyObjFunc<T>>{

    using Base = ode::ObjectiveSolver<ode::RK45, T, 4, ode::SolverPolicy::Static, OdeSystem<T>, MySolver<T>, MyObjFunc<T>>;

public:

    MySolver(T xpoinc, T c1, T c2, T c3, T t0, const T* q0, size_t nsys, T rtol, T atol, T min_step=0, T max_step=ode::inf<T>(), T stepsize=0, int dir=1) : Base(
        ode::ObjFunData<T, MyObjFunc<T>>{.func=MyObjFunc<T>{.point=xpoinc},
                                         .ftol=0, .dir=1}, OdeSystem<T>{.c1=c1, .c2=c2, .c3=c3}, t0, q0, 4, rtol, atol, min_step, max_step, stepsize, dir, std::vector<T>{}){}

};




template <typename Type, typename Type_Arr, typename Type_Vec> class Parent_Poincare
{   
    using Scalar = std::conditional_t<std::is_same_v<Type, MP_REAL>, lazy::LazyType<mpfr::mpreal>, Type>;

  public:
    void print_results(const Type_Arr& results)
    {
        // A function used to print the orbits found by the algorithm in a
        // presentable manner.

        if (results.size() == 0)
        {
            (*this->output) << "/------------------- " << std::endl;
            (*this->output) << "|- No Orbits Found - " << std::endl;
            (*this->output) << "\\-------------------" << std::endl;
        }
        else
        {
            int orbit = 1;

            (*this->output) << std::endl;
            (*this->output) << "/---------------\\" << std::endl;

            for (int i = 0; i < results.cols();)
            {
                (*this->output) << "|--- Orbit " << orbit++ << " ---|" << std::endl;

                for (int j = 0; j < this->pc.p; i++, j++)
                {
                    (*this->output) << "|~ " << results.col(i).transpose() << std::endl;
                }

                (*this->output) << "|--------------";

                for (int j = 0; j < int(std::to_string(orbit).length()); j++)
                {
                    (*this->output) << "-";
                }

                (*this->output) << "|" << std::endl;
            }

            (*this->output) << "\\---------------/" << std::endl;
        }
    }

    Type_Arr calculate(const Type_Vec& particle, Type_Vec &q, bool first_only)
    {
        q[0] = particle(0);
        q[1] = this->pc.xpoin;
        q[2] = particle(1);
        q[3] = 0.0;

        Type_Arr energy_error(particle.rows() + 1, 1);

        for (int i = 0; i < energy_error.rows(); i++)
        {
            energy_error(i, 0) = std::numeric_limits<Type>::infinity();
        }

        if (!(this->check_energy(q.data())))
        {
            // If the energy is out-of-bounds, return error (infinity).

            return energy_error;
        }

        // Type t0 = 0;
        // Type atol = 0;
        // Type rtol = Type(1e-4) * this->pc.err_goal;
        // Type min_step = 0;
        constexpr size_t N = 4;

        std::array<Scalar, N> q0;
        for (size_t i=0; i<N; i++){
            q0[i] = q[i];
        }

        MySolver<Scalar> new_solver(this->pc.xpoin, pc.c1, pc.c2, pc.c3, Scalar(0), q0.data(), 4, Type(1e-4) * pc.err_goal, Type(0), Scalar(0), ode::inf<Type>(), pc.dt);
        // new_solver.set_ics(0, q0.data(), this->pc.dt);

        ode::Array2D<Type, 0, N> ode_results(this->pc.p, N);

        auto failed = [&](const Scalar& t, const Scalar* q) NDSPAN_LAMBDA_INLINE {
            const Scalar& x = q[0];
            const Scalar& y = q[1];
            const Scalar& px = q[2];
            const Scalar& py = q[3];

            Scalar E = (px * px + py * py) / 2 + this->V(x, y);

            Type error = (abs(E - this->pc.ene) / this->pc.ene);
            return error > 1e-6 || t > 1000; // Stop if energy is close enough and we haven't exceeded a reasonable time limit
        };

        int ev_counts = 0;

        while (ev_counts < this->pc.p && new_solver.advance())
        {
            // Store the located poincare section
            if (new_solver.is_at_objective()){
                for (size_t i = 0; i < N; i++){
                    ode_results(ev_counts, i) = new_solver.vector()[i];
                }
                ev_counts++;
            } else if (failed(new_solver.t(), new_solver.vector().data())){
                break; // Stop if energy is close enough and we haven't exceeded a reasonable time limit
            }
        }


        if (ev_counts < this->pc.p)
        {
            return energy_error;
        }
        else if (first_only)
        {
            // If we need only the first poincare section,
            // return the x,px,py values of that section.

            Type_Arr poinc_sections(particle.rows() + 1, 1);

            size_t Nres = this->pc.p;
            for (int j = 0; j <= particle.rows() + 1; j++)
            {

                if (j == 0)
                {
                    poinc_sections(j, 0) = ode_results(Nres - 1, j);
                }
                else if (j > 1)
                {
                    poinc_sections(j - 1, 0) = ode_results(Nres - 1, j);
                }
            }

            return poinc_sections;
        }
        else
        {
            // If we need all the poincare sections found (for p > 1),
            // return all of them.

            Type_Arr poinc_sections(particle.rows(), this->pc.p);

            for (int j = 0; j < particle.rows(); j++)
            {
                for (int k = 0; k < this->pc.p; k++)
                {
                    if (j == 0)
                    {
                        poinc_sections(j, k) = ode_results(k, j);
                    }
                    else if (j >= 1)
                    {
                        poinc_sections(j, k) = ode_results(k, j + 1);
                    }
                }
            }

            return poinc_sections;
        }
    }

  protected:

    Parent_Poincare(const poinc_params<Type> &params, std::ostream *out_stream) : output(out_stream), pc(params){}

    std::ostream *output;
    poinc_params<Type> pc;

#ifdef __MAC__
    dispatch_semaphore_t sem_lock;
#else
    std::mutex mutex_lock;
#endif

    ~Parent_Poincare()
    {
#ifdef __MAC__
        if (this->sem_lock)
        {
            dispatch_release(this->sem_lock);
            this->sem_lock = nullptr;
        }
#endif
    }

    bool check_energy(Type *q)
    {
        Type x = q[0];
        Type y = q[1];
        Type px = q[2];
        // Type py  = q[3];
        Type YY = 2 * (this->pc.ene - this->V(x, y));
        Type YY1 = YY - px*px;

        if (YY1 > 0)
        {
            q[3] = sqrt(YY - px*px);

            return true;
        }
#ifdef __MAC__
        else
        {
            dispatch_semaphore_wait(this->sem_lock, DISPATCH_TIME_FOREVER);

            if (YY1 < 0)
            {
                (*this->output) << "~> Error: The initial energy is not constant" << std::endl;
            }
            else
            {
                (*this->output) << "~> Error: The initial energy is not constant (=0)" << std::endl;
            }

            dispatch_semaphore_signal(this->sem_lock);

            return false;
        }
#else
        else
        {
            std::lock_guard<std::mutex> lock(this->mutex_lock);

            if (YY1 < 0)
            {
                (*this->output) << "~> Error: The initial energy is not constant" << std::endl;
            }
            else
            {
                (*this->output) << "~> Error: The initial energy is not constant (=0)" << std::endl;
            }

            return false;
        }
#endif
    }

    Type V(Type x, Type y) const
    {
        return this->pc.c1 * (x * x + y * y) + this->pc.c2 * y -
               this->pc.c3 * (x * x * x * x + y * y * y * y - 6 * x * x * y * y);
    }

    Type energy(const Type *q) const
    {
        Type x = q[0];
        Type y = q[1];
        Type px = q[2];
        Type py = q[3];

        return ((px * px + py * py) / 2) + this->V(x, y);
    }
};

#endif