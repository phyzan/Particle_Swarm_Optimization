#ifndef __LCL__
#define __LCL__

#include "local_classic_internal.hpp"

class Local_Classic
{
  public:
    Local_Classic(pso_params<MP_REAL, MPFR_ARR> &p, enhanced_params<MP_REAL> &ep, double swap_point, int precision = 64,
                  bool constriction = false, std::ostream *output = &(std::cout));

    ~Local_Classic();

    bool is_initialized() const;
    bool parameter_check();

    void print_params();
    void set_Obj_F(const std::string &objective_func, const obj_params<MP_REAL> &params);

    MPFR_ARR fit();

  private:
    bool initialized;
    double swap_point;

    Local_Classic_Internal<MP_REAL, MPFR_ARR, MPFR_VEC, MPFR_EMPTY> *lclc_in_mpreal;
    Local_Classic_Internal<double, DOUBLE_ARR, DOUBLE_VEC, DOUBLE_EMPTY> *lclc_in_double;
};

inline bool Local_Classic::is_initialized() const
{
    return this->initialized;
};

inline bool Local_Classic::parameter_check()
{
    return this->lclc_in_mpreal->parameter_check();
};

inline void Local_Classic::print_params()
{
    this->lclc_in_mpreal->print_params();
}

Local_Classic *Local_Classic_Init(pso_params<MP_REAL, MPFR_ARR> &p, enhanced_params<MP_REAL> &ep, double swap_point,
                                  int precision = 64, bool constriction = false, std::ostream *output = &(std::cout));

#endif