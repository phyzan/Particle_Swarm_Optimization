#ifndef __G_H__
#define __G_H__

#include "../../../functions/functions.hpp"
#include "../../../local_definitions.hpp"

template <typename Type, typename Type_Vec> class H
{
  public:
    H(int d, int w)
    {
        this->d = d;
        this->w = w;
        this->b = uniform_real(0, w);
        this->v = Type_Vec::NullaryExpr(d, [&]() { return Type(gauss_distribution(0, 1)); });
    }

    uint64_t calculate(Type_Vec particle)
    {
        return uint64_t(lmath::floor((particle.array().cwiseProduct(this->v.array()).sum() + this->b) / this->w));
    }

  private:
    int d;      // The dimension of the particle.
    int w;      // A "window" size.
    Type b;     // A real number sampled from [0,w).
    Type_Vec v; // A random vector initialized with the Gaussian distribution.
};

template <typename Type, typename Type_Vec> class G
{
  public:
    G(int d, int w, int k)
    {
        this->k = k;
        this->M = pow(2, 32) - 5;
        this->h_ = new H<Type, Type_Vec> *[k];
        this->r_ = new uint64_t[k];

        for (int i = 0; i < this->k; i++)
        {
            this->h_[i] = new H<Type, Type_Vec>(d, w);
            this->r_[i] = uniform_int(0, this->M);
        }
    }

    ~G()
    {
        for (int i = 0; i < this->k; i++)
        {
            delete this->h_[i];

            this->h_[i] = nullptr;
        }

        delete[] this->h_;
        delete[] this->r_;

        this->h_ = nullptr;
        this->r_ = nullptr;
    }

    uint64_t calculate_id(Type_Vec particle)
    {
        uint64_t result = 0;

        for (int i = 0; i < this->k; i++)
        {
            result += this->r_[i] * this->h_[i]->calculate(particle);
        }

        return result % this->M;
    }

    uint64_t calculate_g(Type_Vec particle, int table_size)
    {
        return this->calculate_id(particle) % table_size;
    }

  private:
    int k;        // Number of h() functions used for the G() function.
    uint64_t M;   // A large integer number.
    uint64_t *r_; // Random integer numbers.

    H<Type, Type_Vec> **h_;
};

#endif