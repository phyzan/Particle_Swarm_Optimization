#ifndef _SPNC_
#define _SPNC_

#include "parent_poincare.hpp"

template <typename Type, typename Type_Arr, typename Type_Vec>
class Simple_Poincare : public Parent_Poincare<Type, Type_Arr, Type_Vec>
{
  public:
    Simple_Poincare(const poinc_params<Type> &pc, std::ostream *output = &(std::cout))
    {
        this->pc = pc;
        this->output = output;

#ifdef __MAC__
        this->sem_lock = dispatch_semaphore_create(1);
#endif
    }

    Type_Vec objective_function(Type_Arr population)
    {
        Type_Vec fpopulation(population.cols());
#ifdef __MAC__
        if (this->pc.threads > 1)
        { // Multithreading for MAC OS.
            dispatch_group_t group = dispatch_group_create();
            dispatch_queue_t queue = dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0);

            for (int i = 0; i < population.cols(); i++)
            {
                InstanceParams *instance = new InstanceParams{i, this, &fpopulation, population};

                dispatch_group_async_f(group, queue, instance, Simple_Poincare::objf_calculation);
            }

            dispatch_group_wait(group, DISPATCH_TIME_FOREVER);
        }
#else
        if (this->pc.threads > 1)
        { // Multithreading for Linux.
            omp_set_num_threads(this->pc.threads);
            Eigen::setNbThreads(this->pc.threads);

#pragma omp parallel for
            for (int i = 0; i < population.cols(); i++)
            {
                InstanceParams *instance = new InstanceParams{i, this, &fpopulation, population};

                Simple_Poincare::objf_calculation(instance);
            }
        }
#endif
        else
        { // No multithreading (user declared pc.threads <= 1).
            for (int i = 0; i < population.cols(); i++)
            {
                InstanceParams *instance = new InstanceParams{i, this, &fpopulation, population};

                Simple_Poincare::objf_calculation(instance);
            }
        }

        return fpopulation;
    }

  private:
    struct InstanceParams
    {
        int i; // The number of particle that is being processed.
        Simple_Poincare *instance;
        Type_Vec *fpopulation;
        const Type_Arr population;
    };

    static void objf_calculation(void *instance)
    {
        Type_Vec q(4);
        InstanceParams *l_instance = static_cast<InstanceParams *>(instance);

        if (lmath::isinf(l_instance->population(0, l_instance->i)))
        {
            // If the current particle has out-of-bounds energy (i.e. =inf), do not process it and return.

            delete l_instance;
            l_instance = nullptr;
            return;
        }

        // Calculate the new position of the particle.
        Type_Arr poinc_sections = l_instance->instance->calculate((l_instance->population.col(l_instance->i)), q, true);

        Type x = poinc_sections(0, poinc_sections.cols() - 1);
        Type px = poinc_sections(1, poinc_sections.cols() - 1);
        Type py = poinc_sections(2, poinc_sections.cols() - 1);

        // Calculate the difference of the new and the old position of the particle using
        // the Euclidean Distance.
        if (py > 0 && x != lmath::get_infinity<Type>())
        {
            (*l_instance->fpopulation)(l_instance->i) = lmath::sqrt(lmath::pow(q[0] - x, 2) + lmath::pow(q[2] - px, 2));
        }
        else
        { // If py <= 0, return infinity.
            (*l_instance->fpopulation)(l_instance->i) = lmath::get_infinity<Type>();
        }

        delete l_instance;
        l_instance = nullptr;
    }
};

#endif