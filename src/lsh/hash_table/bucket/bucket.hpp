#ifndef __BCKT__
#define __BCKT__

#include "../../../local_definitions.hpp"

template <typename Type_Vec> class item
{
  public:
    item(int index, Type_Vec particle)
    {
        this->next = nullptr;
        this->index = index;
        this->particle = particle;
    }

    ~item()
    {
        this->next = nullptr;
    }

    inline int get_index()
    {
        return this->index;
    }

    inline item *get_next()
    {
        return this->next;
    }

    inline void set_next(item *next)
    {
        this->next = next;
    }

    inline Type_Vec get_particle()
    {
        return this->particle;
    }

  private:
    int index; // The index of the particle in the population array of PSO.
    Type_Vec particle;

    item *next;
};

template <typename Type_Vec> class bucket
{
  public:
    bucket()
    {
        this->size = 0;
        this->item_list = nullptr;
    }

    ~bucket()
    {
        while (this->item_list)
        {
            item<Type_Vec> *temp = this->item_list;
            this->item_list = this->item_list->get_next();

            delete temp;
            temp = nullptr;
        }
    }

    void add_particle(int index, Type_Vec particle)
    {
        item<Type_Vec> *new_item = new item<Type_Vec>(index, particle);

        new_item->set_next(this->item_list); // Add the new item at the head of the list.
        this->item_list = new_item;
        this->size++; // Increase the number of items in the bucket.
    }

    inline item<Type_Vec> *get_items()
    {
        return this->item_list;
    }

  private:
    int size; // The number of items in the hash_bucket.
    item<Type_Vec> *item_list;
};

#endif