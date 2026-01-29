#ifndef __LSH__
#define __LSH__

#include "hash_table/hash_table.hpp"
#include "priority_list/priority_list.hpp"

template <typename Type, typename Type_Arr, typename Type_Vec> class LSH
{
  public:
    LSH(int k, int w, int L, Type err_goal, Type effective_radius, Type_Arr dim_limits, Type_Arr population)
    {
        this->L = L;
        this->err_goal = 1e-4 * err_goal;
        this->N = population.cols() / 4;
        this->dim_limits = dim_limits;
        this->effective_radius = effective_radius;

        this->HT = new hash_table<Type, Type_Arr, Type_Vec> *[this->L];

        for (int i = 0; i < this->L; i++)
        {
            G<Type, Type_Vec> *g = new G<Type, Type_Vec>(population.rows(), w, k);
            this->HT[i] = new hash_table<Type, Type_Arr, Type_Vec>(g, population);
        }
    }

    ~LSH()
    {
        for (int i = 0; i < this->L; i++)
        {
            delete this->HT[i];
            this->HT[i] = nullptr;
        }
        delete[] this->HT;
        this->HT = nullptr;
    }

    priority_list<Type, Type_Vec> *find(Type_Vec particle)
    {
        priority_list<Type, Type_Vec> *neighbours = new priority_list<Type, Type_Vec>(this->err_goal, particle);
        // Create a priority list for the neighbours. A priority list is needed
        // because if/when the list reaches N+1 items, the farthest item should be
        // removed. With a priority list, that item will always be the last item,
        // making its removal time-efficient.

        for (int i = 0; i < this->L; i++)
        { // For each hash table,
            item<Type_Vec> *bucket_items =
                this->HT[i]->get_bucket_items(particle); // get the items from the bucket that
                                                         // the particle would be a part of,
            uint64_t particle_id =
                this->HT[i]->get_hash_function()->calculate_id(particle); // calculate the particle id for the
                                                                          // specific hash table,
            while (bucket_items)
            {
                int bucket_particle_idx = bucket_items->get_index();
                Type_Vec bucket_particle(bucket_items->get_particle());

                uint64_t item_id = this->HT[i]->get_hash_function()->calculate_id(
                    bucket_particle); // and calculate the current item's hash id.

                if (particle_id == item_id)
                { // If the IDs match,
                    if (this->effective_radius == 1)
                    { // if the radius is 1,
                        if (neighbours->particle_exists(bucket_particle) == false)
                        { // and the item is not already in the neighbours list,
                            neighbours->add_particle(bucket_particle_idx, bucket_particle); // add it to the list.
                        }
                    }
                    else
                    { // If effective_radius < 1, then calculate the new
                      // search space.
                        Type_Vec max_dim_distance = Type_Vec(this->dim_limits.col(1) - this->dim_limits.col(0));
                        Type_Vec max_particle_radius = particle + ((this->effective_radius / 2) * max_dim_distance);
                        Type_Vec min_particle_radius = particle - ((this->effective_radius / 2) * max_dim_distance);

                        if ((bucket_particle.array() >= min_particle_radius.array() &&
                             bucket_particle.array() < max_particle_radius.array())
                                .all())
                        { // If the item is in the effective radius
                            if (neighbours->particle_exists(bucket_particle) == false)
                            { // and it is not already in the neighbours list,
                                neighbours->add_particle(bucket_particle_idx, bucket_particle); // add it to the list.
                            }
                        }
                    }

                    if (neighbours->get_size() > this->N)
                    { // If the neighbours list has more than N items,
                      // remove the farthest (=last) one.
                        neighbours->remove_particle();
                    }
                }

                bucket_items = bucket_items->get_next();
            }
        }

        return neighbours;
    }

    Eigen::Vector<int, Eigen::Dynamic> find_idx(Type_Vec particle)
    {
        priority_list<Type, Type_Vec> *neighbours =
            this->find(particle); // Find the N closest neighbours of the particle,

        Eigen::Vector<int, Eigen::Dynamic> idx(neighbours->get_size()); // create an integer array and

        int i = 0;
        priority_node<Type, Type_Vec> *temp = neighbours->get_items();

        while (temp)
        {
            idx(i) = temp->get_index(); // save all the indices of the neighbours.

            i++;
            temp = temp->get_next();
        }

        delete neighbours;
        neighbours = nullptr;

        return idx;
    }

  private:
    int L; // Number of hash tables.
    int N; // Number of nearest neighbours the algorithm returns.
    Type err_goal;
    Type effective_radius; // With the query particle as the center,
                           // the algorithm creates a hypercube that
                           // is a percentage of the size of the initial space
                           // limits (dim_limits). If effective_radius = 1,
                           // the initial limits (dim_limits) are used and
                           // the query particle is not treated as the center.
    Type_Arr dim_limits;   // The initial limits of the search space.

    hash_table<Type, Type_Arr, Type_Vec> **HT;
};

#endif