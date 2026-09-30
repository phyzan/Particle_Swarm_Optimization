#ifndef PSO_INTEGRATE_SYSTEM_HPP
#define PSO_INTEGRATE_SYSTEM_HPP


#include <src/input/parameters.hpp>
#include <odecraft/odecraft.hpp>

using xdiff::power;

/// @brief Hamiltonian of a caldera-like potential energy surface.
template<typename T>
struct CalderaODE{

    T c1, c2, c3;

    /// @brief Potential energy of the caldera surface at (x, y).
    XDIFF_FORCEINLINE decltype(auto) V(const auto& x, const auto& y) const{
        return c1*(power<2>(x) + power<2>(y)) + c2*y - c3*(power<4>(x) + power<4>(y) - 6*power<2>(x)*power<2>(y));
        // return c1*(x*x + y*y) + c2*y
            //  - c3*(x*x*x*x + y*y*y*y - 6*x*x*y*y);
    }

    /// @brief Hamilton's equations: writes dq/dt for the state q.
    XDIFF_FORCEINLINE void Rhs(T* out, const T& t, const T* q) const {
        // auto qq = xdiff::as_seedvector<2, 1, xdiff::Layout::Nested>(q);
        
        // const auto &x = qq[0];
        // const auto &y = qq[1];
        // const T &px = q[2];
        // const T &py = q[3];

        // auto v = this->V(x, y);
        // out[0] = px;
        // out[1] = py;
        // out[2] = -v.get_diff_wrt(0);
        // out[3] = -v.get_diff_wrt(1);

        const T& x = q[0];
        const T& y = q[1];
        const T& px = q[2];
        const T& py = q[3];
        auto vx = 2 * c1 * x - 4 * c3 * x * x * x + 12 * c3 * x * y * y;
        auto vy = 2 * c1 * y + c2 - 4 * c3 * y * y * y + 12 * c3 * y * x * x;

        out[0] = px;
        out[1] = py;
        out[2] = -vx;
        out[3] = -vy;
    }
};

/// @brief Builds the CalderaODE the parameters describe.
template<typename T>
inline CalderaODE<T> make_system(const input::Parameters& pin)
{
    return CalderaODE<T>{
        .c1 = T{pin.objective.poincare.c1},
        .c2 = T{pin.objective.poincare.c2},
        .c3 = T{pin.objective.poincare.c3},
    };
}

#endif // PSO_INTEGRATE_SYSTEM_HPP