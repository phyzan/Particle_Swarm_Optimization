#ifndef __PRND__
#define __PRND__

#include "../../../local_definitions.hpp"

template <typename Type, typename Type_Vec> class priority_node
{
  public:
    priority_node(int index, Type_Vec particle, Type dist)
    {
        // Initialize the pointers to their initial values.

        this->dist = dist;
        this->next = nullptr;
        this->prev = nullptr;
        this->index = index;
        this->particle = particle;
    }

    ~priority_node()
    {
        this->next = nullptr;
        this->prev = nullptr;
    }

    inline int get_index()
    {
        return this->index;
    }

    inline Type get_dist()
    {
        return this->dist;
    }

    inline Type_Vec get_particle()
    {
        return this->particle;
    }

    inline priority_node *get_next()
    {
        return this->next;
    }

    inline priority_node *get_prev()
    {
        return this->prev;
    }

    inline void set_next(priority_node *next)
    {
        this->next = next;
    }

    inline void set_prev(priority_node *prev)
    {
        this->prev = prev;
    }

  private:
    int index;         // The index of the particle in the population array of PSO.
                       // It is used by the LSH algorithm to return indices to the
                       // population list of the PSO algorithm instead of particles
                       // in the population array, to conserve time.
    Type dist;         // The distance of the node particle from the base particle.
    Type_Vec particle; // The particle of the node.

    priority_node *next;
    priority_node *prev;
};

#endif