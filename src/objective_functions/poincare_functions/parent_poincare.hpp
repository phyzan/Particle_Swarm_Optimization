#ifndef __PPNC__
#define __PPNC__

#include "../../odepack/include/odepack/ode.hpp"
#include <thread>

#ifdef __MAC__
#include <dispatch/dispatch.h>
#else
#include <mutex>
#include <omp.h>
#endif

#include "../../local_definitions.hpp"

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

template <typename Type> void F(Type *res, const Type &t, const Type *q, const Type *args, const void *aux)
{
    // The function F that calculates the new position.
    const auto *pc = reinterpret_cast<const poinc_params<Type> *>(aux);

    const Type &x = q[0];
    const Type &y = q[1];
    const Type &px = q[2];
    const Type &py = q[3];

    Type vx = 2 * pc->c1 * x - 4 * pc->c3 * x * x * x + 12 * pc->c3 * x * y * y;
    Type vy = 2 * pc->c1 * y + pc->c2 - 4 * pc->c3 * y * y * y + 12 * pc->c3 * y * x * x;

    res[0] = px;
    res[1] = py;
    res[2] = -vx;
    res[3] = -vy;
};

template <typename Type> void Jac(Type *J, const Type &t, const Type *q, const Type *args, const void *aux)
{
    const auto *pc = reinterpret_cast<const poinc_params<Type> *>(aux);

    const Type &x = q[0];
    const Type &y = q[1];

    // Precompute terms
    Type dvx_dx = 2 * pc->c1 - 12 * pc->c3 * x * x + 12 * pc->c3 * y * y;
    Type dvx_dy = 24 * pc->c3 * x * y;

    Type dvy_dx = 24 * pc->c3 * x * y;
    Type dvy_dy = 2 * pc->c1 - 12 * pc->c3 * y * y + 12 * pc->c3 * x * x;

    // Fill Jacobian in row-major order:
    // J[i*4 + j] = dF_i / dq_j

    // F0 = px
    J[0 * 4 + 0] = Type(0);
    J[0 * 4 + 1] = Type(0);
    J[0 * 4 + 2] = Type(1);
    J[0 * 4 + 3] = Type(0);

    // F1 = py
    J[1 * 4 + 0] = Type(0);
    J[1 * 4 + 1] = Type(0);
    J[1 * 4 + 2] = Type(0);
    J[1 * 4 + 3] = Type(1);

    // F2 = -vx
    J[2 * 4 + 0] = -dvx_dx;
    J[2 * 4 + 1] = -dvx_dy;
    J[2 * 4 + 2] = Type(0);
    J[2 * 4 + 3] = Type(0);

    // F3 = -vy
    J[3 * 4 + 0] = -dvy_dx;
    J[3 * 4 + 1] = -dvy_dy;
    J[3 * 4 + 2] = Type(0);
    J[3 * 4 + 3] = Type(0);
}

template <typename T> T _event(const T &t, const T *q, const T *args, const void *aux)
{
    // The event function, it is true when it is close to 0 (meaning y ~= x_poin).
    const auto *pc = reinterpret_cast<const poinc_params<T> *>(aux);

    return q[1] - pc->xpoin;
}

template <typename T> T _stop_event(const T &t, const T *q, const T *args, const void *aux)
{
    // A function that stops the ODE calculations if the energy goes out-of-bounds.
    const auto *p = reinterpret_cast<const Parent_Poincare<T, Eigen::Array<T, -1, -1>, Eigen::Vector<T, -1>> *>(aux);

    T E = p->energy(q);

    return (lmath::abs(E - p->pc.ene) / p->pc.ene) - T(1e-6);
}

template <typename Type> class StopEvent : public PreciseEvent<Type, static_cast<size_t>(4)>
{ // A custom ODE stop event.
  public:
    StopEvent(std::string name, ObjFun<Type> f, const void *pp)
        : PreciseEvent<Type, static_cast<size_t>(4)>(name, f, 1, nullptr, false, 1e-10, pp)
    {
    }

    bool is_leathal() const override
    {
        return true;
    };
    bool is_stop_event() const override
    {
        return true;
    };

    StopEvent *clone() const override
    {
        return new StopEvent(*this);
    };
};

template <typename Type, typename Type_Arr, typename Type_Vec> class Parent_Poincare
{
    friend Type _stop_event<Type>(const Type &t, const Type *q, const Type *args, const void *aux);

  public:
    void print_results(Type_Arr results)
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

    Type_Arr calculate(Type_Vec particle, Type_Vec &q, bool first_only)
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

        if (this->check_energy(q.data()) == false)
        {
            // If the energy is out-of-bounds, return error (infinity).

            return energy_error;
        }

        Type t0 = 0;
        Type atol = 0;
        Type rtol = Type(1e-4) * this->pc.err_goal;
        Type min_step = 0;
        constexpr size_t N = 4;

        StopEvent<Type> stop_events("poincare_stop", _stop_event<Type>, this);
        PreciseEvent<Type, N> events("poincare_sect", _event<Type>, 1, nullptr, false, 1e-20, &(this->pc));

        Array1D<Type, N> q0(q.data());

        RK45<Type, N> solver({F, Jac, &(this->pc)}, t0, q0, rtol, atol, min_step, inf<Type>(), this->pc.dt, 1, {},
                             {&events, &stop_events});

        Array2D<Type, 0, N> ode_results(this->pc.p, N);

        int ev_counts = 0;
        Type tmax = 1000;

        while (ev_counts < this->pc.p && solver.t() < tmax && !solver.is_dead() && solver.advance_to_event())
        {
            for (size_t ev : solver.event_col())
            {
                if (ev == 0)
                {
                    copy_array(ode_results.data() + N * ev_counts, solver.q().data(), N);
                    ev_counts++;
                }
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
        Type YY1 = YY - lmath::pow(px, 2);

        if (YY1 > 0)
        {
            q[3] = lmath::sqrt(YY - lmath::pow(px, 2));

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
        return this->pc.c1 * (lmath::pow(x, 2) + lmath::pow(y, 2)) + this->pc.c2 * y -
               this->pc.c3 * (lmath::pow(x, 4) + lmath::pow(y, 4) - 6 * lmath::pow(x, 2) * lmath::pow(y, 2));
    }

    Type energy(const Type *q) const
    {
        Type x = q[0];
        Type y = q[1];
        Type px = q[2];
        Type py = q[3];

        return ((lmath::pow(px, 2) + lmath::pow(py, 2)) / 2) + this->V(x, y);
    }
};

#endif