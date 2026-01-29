#ifndef __GCL__
#define __GCL__

#include "global_classic_internal.hpp"

class Global_Classic
{
  public:
    Global_Classic(pso_params<MP_REAL, MPFR_ARR> &p, double swap_point, int precision = 64,
                   std::ostream *output = &(std::cout));
    ~Global_Classic();

    bool is_initialized() const;
    bool parameter_check();

    void print_params();
    void set_Obj_F(const std::string &objective_func, const obj_params<MP_REAL> &params);

    MPFR_ARR fit();

  private:
    bool initialized;
    double swap_point;

    Global_Classic_Internal<MP_REAL, MPFR_ARR, MPFR_VEC, MPFR_EMPTY> *gclc_in_mpreal;
    Global_Classic_Internal<double, DOUBLE_ARR, DOUBLE_VEC, DOUBLE_EMPTY> *gclc_in_double;
};

inline bool Global_Classic::is_initialized() const
{
    return this->initialized;
};

inline bool Global_Classic::parameter_check()
{
    return this->gclc_in_mpreal->parameter_check();
};

inline void Global_Classic::print_params()
{
    this->gclc_in_mpreal->print_params();
}

Global_Classic *Global_Classic_Init(pso_params<MP_REAL, MPFR_ARR> &p, double swap_point, int precision = 64,
                                    std::ostream *output = &(std::cout));

#endif