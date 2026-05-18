#ifndef _PNC_
#define _PNC_

#include "parent_poincare.hpp"

template <typename Type, typename Type_Arr, typename Type_Vec>
class Poincare : public Parent_Poincare<Type, Type_Arr, Type_Vec>
{
  public:
    Poincare(const poinc_params<Type> &pc, std::ostream *output = &(std::cout)) : Parent_Poincare<Type, Type_Arr, Type_Vec>(pc, output)
    {
#ifdef __MAC__
        this->sem_lock = dispatch_semaphore_create(1);
#endif
    }

    Type_Vec objective_function(const Type_Arr& population)
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

                dispatch_group_async_f(group, queue, instance, Poincare::objf_calculation);
            }

            dispatch_group_wait(group, DISPATCH_TIME_FOREVER);
            dispatch_release(group);
        }
#else
        if (this->pc.threads > 1)
        { // Multithreading for Linux.
            omp_set_num_threads(this->pc.threads);

#pragma omp parallel for schedule(dynamic)
            for (int i = 0; i < population.cols(); i++)
            {
                // Stack allocation instead of heap - no new/delete overhead
                objf_calculation_inline(i, this, fpopulation, population);
            }
        }
#endif
        else
        { // No multithreading (user declared pc.threads <= 1).
            for (int i = 0; i < population.cols(); i++)
            {
                objf_calculation_inline(i, this, fpopulation, population);
            }
        }

        return fpopulation;
    }

  private:
    struct InstanceParams
    {
        int i; // The number of particle that is being processed.
        Poincare *instance;
        Type_Vec *fpopulation;
        const Type_Arr population;
    };

    // inline calculation without heap allocation - used by OpenMP and single-threaded paths
    static void objf_calculation_inline(int i, Poincare *instance, Type_Vec &fpopulation, const Type_Arr &population)
    {
        Type_Vec q(4);

        if (isinf(population(0, i)))
        {
            // If the current particle has out-of-bounds energy (i.e. =inf), do not process it.
            return;
        }

        // Calculate the new position of the particle.
        Type_Arr poinc_sections = instance->calculate((population.col(i)), q, true);

        Type x = poinc_sections(0, poinc_sections.cols() - 1);
        Type px = poinc_sections(1, poinc_sections.cols() - 1);
        Type py = poinc_sections(2, poinc_sections.cols() - 1);

        Type magnitude_p = (x*x + px*px + py*py);

        // Calculate the difference of the new and the old position of the particle using the
        // function described in the publication.
        if (py > 0 && x != ode::inf<Type>())
        {
            fpopulation(i) = sqrt((x * x * (q[0] - x) * (q[0] - x)) +
                                  (px * px * (q[2] - px) * (q[2] - px)) +
                                  (py * py * (q[3] - py) * (q[3] - py))) /
                             magnitude_p;
        }
        else
        { // If py <= 0, return infinity.
            fpopulation(i) = ode::inf<Type>();
        }

        mpfr_free_cache2(MPFR_FREE_LOCAL_CACHE);
    }

    // Heap-based version for Mac dispatch (requires void* parameter)
    static void objf_calculation(void *instance)
    {
        Type_Vec q(4);
        InstanceParams *l_instance = static_cast<InstanceParams *>(instance);

        if (isinf(l_instance->population(0, l_instance->i)))
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

        Type magnitude_p = (x * x + px * px + py * py);

        // Calculate the difference of the new and the old position of the particle using the
        // function described in the publication.
        if (py > 0 && x != ode::inf<Type>())
        {
            (*l_instance->fpopulation)(l_instance->i) = sqrt((x * x * (q[0] - x) * (q[0] - x)) +
                                                                    (px * px * (q[2] - px) * (q[2] - px)) +
                                                                    (py * py * (q[3] - py) * (q[3] - py))) /
                                                        magnitude_p;
        }
        else
        { // If py <= 0, return infinity.
            (*l_instance->fpopulation)(l_instance->i) = ode::inf<Type>();
        }

        delete l_instance;
        l_instance = nullptr;

        mpfr_free_cache2(MPFR_FREE_LOCAL_CACHE);
    }
};

#endif