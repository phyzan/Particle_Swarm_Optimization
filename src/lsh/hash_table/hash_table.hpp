#ifndef __HASH__
#define __HASH__

#include "bucket/bucket.hpp"
#include "g_h/g_h.hpp"

template <typename Type, typename Type_Arr, typename Type_Vec> class hash_table
{
  public:
    hash_table(G<Type, Type_Vec> *g, Type_Arr population)
    {
        int n = population.cols();

        this->g = g;
        this->bucket_num = ceil(n / 8); // Number of buckets.
        this->bucket_ = new bucket<Type_Vec> *[this->bucket_num];

        for (int i = 0; i < this->bucket_num; i++)
        {
            this->bucket_[i] = new bucket<Type_Vec>;
        }

        for (int i = 0; i < n; i++)
        {
            // For each particle calculate its hash and place it
            // in the respective bucket.

            int g_id = this->g->calculate_g(population.col(i), this->bucket_num);

            this->bucket_[g_id]->add_particle(i, population.col(i));
        }
    }

    ~hash_table()
    {
        for (int i = 0; i < this->bucket_num; i++)
        {
            delete this->bucket_[i];
            this->bucket_[i] = nullptr;
        }
        delete[] this->bucket_;
        this->bucket_ = nullptr;

        delete this->g;
        this->g = nullptr;
    }

    item<Type_Vec> *get_bucket_items(Type_Vec particle)
    {
        // Return the items stored in a specific bucket.

        int g_id = this->g->calculate_g(particle, this->bucket_num);

        return this->bucket_[g_id]->get_items();
    }

    inline G<Type, Type_Vec> *get_hash_function()
    {
        return this->g;
    };

  private:
    G<Type, Type_Vec> *g;       // The hash function.
    bucket<Type_Vec> **bucket_; // An array of the hash buckets.
    int bucket_num;             // Number of buckets in the hash table.
};

#endif