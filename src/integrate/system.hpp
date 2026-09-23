#ifndef PSO_INTEGRATE_SYSTEM_HPP
#define PSO_INTEGRATE_SYSTEM_HPP

#include <src/input/parameters.hpp>

#include <odecraft/odecraft.hpp>

/// @brief Hamiltonian of a caldera-like potential energy surface.
template<typename T>
struct CalderaODE{

    T c1, c2, c3;

    /// @brief Potential energy of the caldera surface at (x, y).
    inline T V(const T& x, const T& y) const{
        return c1*(x*x + y*y) + c2*y
             - c3*(x*x*x*x + y*y*y*y - 6*x*x*y*y);
    }

    /// @brief Hamilton's equations: writes dq/dt for the state q.
    XDIFF_FORCEINLINE void Rhs(auto* out, const auto& /*t*/, auto q) const{

        const auto &x = q[0];
        const auto &y = q[1];
        const auto &px = q[2];
        const auto &py = q[3];

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
        .c3 = T{pin.objective.poincare.c3}
    };
}

#endif // PSO_INTEGRATE_SYSTEM_HPP
