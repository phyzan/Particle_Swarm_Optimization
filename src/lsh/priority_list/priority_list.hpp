#ifndef __PRIO__
#define __PRIO__

#include "priority_node/priority_node.hpp"

template <typename Type, typename Type_Vec> class priority_list
{
  public:
    priority_list(Type err_goal, Type_Vec base_particle)
    {
        this->size = 0;
        this->err_goal = err_goal;
        this->list_head = nullptr;
        this->list_tail = nullptr;
        this->base_particle = base_particle;
    }

    ~priority_list()
    {
        while (this->list_tail)
        {
            this->remove_particle();
        }

        this->list_head = nullptr;
        this->list_tail = nullptr;
    }

    bool particle_exists(Type_Vec particle)
    {
        // To find if a particle exists in the list, insteal of
        // finding the exact same particle (norm(particle - curr) = 0),
        // find if the norm of the particle with an item in the list
        // is < some_error.
        // If MPreal numbers are not used, this ensures that floating
        // point inaccuraces won't prevent the program from functioning.

        priority_node<Type, Type_Vec> *curr = this->list_head;

        while (curr)
        {
            if ((particle - curr->get_particle()).norm() <= this->err_goal)
            {
                return true;
            }

            curr = curr->get_next();
        }

        return false;
    }

    void remove_particle()
    {
        if (this->list_tail)
        {
            // If there exists a node at the end of the list,
            // remove it.

            if (this->list_tail == this->list_head)
            {
                // If it is the last node in the list,
                // remove it and make the head and tail
                // point to nullptr.

                delete this->list_head;
                this->list_head = nullptr;
                this->list_tail = nullptr;
            }
            else
            {
                priority_node<Type, Type_Vec> *temp = this->list_tail;

                this->list_tail = this->list_tail->get_prev();

                this->list_tail->set_next(nullptr);

                delete temp;
                temp = nullptr;
            }

            this->size--; // Reduce the size of the list by 1.
        }
    }

    void add_particle(int index, Type_Vec particle)
    {
        // To add an item to the list,

        Type new_dist = (this->base_particle - particle).norm();
        // calculate its distance to the base_particle,

        priority_node<Type, Type_Vec> *new_node = new priority_node<Type, Type_Vec>(index, particle, new_dist);
        // create a new list node that contains the index, the particle and the distance,

        if (this->list_head == nullptr)
        {
            // check if the head (and tail) are empty,
            // and if they are, add the new node there.

            this->list_head = new_node;
            this->list_tail = new_node;
        }
        else
        {
            // If they are not,
            priority_node<Type, Type_Vec> *curr = this->list_head;

            while (1)
            {
                // iterate through the list.

                if (new_dist < curr->get_dist())
                {
                    // If the distance of the new node is smaller than that
                    // of the current node, add the new node in front of the
                    // current node.

                    if (this->list_head == curr)
                    {
                        // If the current node is the list_head, make the new
                        // node the list_head.

                        this->list_head = new_node;
                    }
                    else
                    {
                        priority_node<Type, Type_Vec> *prev = curr->get_prev();

                        prev->set_next(new_node);
                        new_node->set_prev(prev);
                    }

                    new_node->set_next(curr);
                    curr->set_prev(new_node);

                    break;
                }

                if (curr->get_next() == nullptr)
                {
                    // If the new node is not closer to the base_particle than
                    // any other node in the list, add it at the end of the list.

                    new_node->set_prev(curr);
                    curr->set_next(new_node);

                    this->list_tail = new_node;

                    break;
                }

                curr = curr->get_next();
            }
        }

        this->size++; // Increase the list size by 1.
    }

    inline int get_size()
    {
        return this->size;
    }

    inline priority_node<Type, Type_Vec> *get_items()
    {
        return this->list_head;
    }

  private:
    int size;                                 // The number of particles in the list.
    Type err_goal;                            // The error allowed for two particles to be considered
                                              // the "same".
    Type_Vec base_particle;                   // The particle that the distances of all the
                                              // other particles are measured from.
    priority_node<Type, Type_Vec> *list_head; // The beginning of the list.
    priority_node<Type, Type_Vec> *list_tail; // The end of the list.
};

#endif