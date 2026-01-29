#ifndef __DEF__
#define __DEF__

#include "deflection_internal.hpp"

class Deflection
{
  public:
    Deflection(pso_params<MP_REAL, MPFR_ARR> &p, enhanced_params<MP_REAL> &ep, double swap_point, MP_REAL &l,
               MP_REAL &constant, int precision = 64, bool constriction = false, bool LLSH = false,
               bool repulsion = false, std::ostream *output = &(std::cout));
    ~Deflection();

    bool is_initialized() const;
    bool parameter_check();

    void print_params();
    void set_Obj_F(const std::string &objective_func, const obj_params<MP_REAL> &params);

    MPFR_ARR fit(int runs = 1);

  private:
    bool initialized;
    double swap_point;

    Deflection_Internal<MP_REAL, MPFR_ARR, MPFR_VEC, MPFR_EMPTY> *defl_in_mpreal;
    Deflection_Internal<double, DOUBLE_ARR, DOUBLE_VEC, DOUBLE_EMPTY> *defl_in_double;
};

inline bool Deflection::is_initialized() const
{
    return this->initialized;
};

inline bool Deflection::parameter_check()
{
    return this->defl_in_mpreal->parameter_check();
};

inline void Deflection::print_params()
{
    this->defl_in_mpreal->print_params();
}

Deflection *Deflection_Init(pso_params<MP_REAL, MPFR_ARR> &p, enhanced_params<MP_REAL> &ep, double swap_point,
                            MP_REAL l, MP_REAL constant, int precision = 64, bool constriction = false,
                            bool LLSH = false, bool repulsion = false, std::ostream *output = &(std::cout));

#endif