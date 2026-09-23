#ifndef PSO_OBJECTIVE_ANALYTIC_HPP
#define PSO_OBJECTIVE_ANALYTIC_HPP

#include <src/objective/objective.hpp>

#include <cmath>
#include <string>

/// @brief Standard analytic test functions, for exercising the swarm itself.
///
/// Evaluation is a handful of flops rather than an ODE integration, so a run
/// finishes in milliseconds. Useful for checking convergence behaviour,
/// topologies and the precision swap without waiting on the physics -- and as
/// the worked example of what a second objective has to provide.
///
/// All three have their global minimum at value 0:
///   sphere      x = 0
///   rastrigin   x = 0
///   rosenbrock  x = 1
template<typename T>
class AnalyticObjective final : public Objective<T>{
public:

    enum class Function : uint8_t{ Sphere, Rastrigin, Rosenbrock };

    /// @brief Selects the function the parameters name.
    explicit AnalyticObjective(const input::Parameters& pin)
        : which(parse(pin.objective.analytic.function)){}

    /// @brief Value at x. Always finite: every point is feasible.
    T evaluate(const T* x) override{

        switch (which){

            case Function::Sphere:{
                T sum = T{0};
                for (size_t i = 0; i < PSO_DIM; i++){ sum += x[i]*x[i]; }
                return sum;
            }

            case Function::Rastrigin:{
                using std::cos;
                T sum = T{10 * PSO_DIM};
                for (size_t i = 0; i < PSO_DIM; i++){
                    sum = sum + x[i]*x[i] - T{10}*cos(2*M_PI*x[i]);
                }
                return sum;
            }

            case Function::Rosenbrock:{
                T sum = T{0};
                for (size_t i = 0; i + 1 < PSO_DIM; i++){
                    const T a = T{1} - x[i];
                    const T b = x[i+1] - x[i]*x[i];
                    sum = sum + a*a + 100*b*b;
                }
                return sum;
            }
        }

        return ode::inf<T>();
    }

    /// @brief The minimiser itself; there is nothing further to reconstruct.
    std::vector<std::array<T, PSO_DIM>> refine(const T* x) override{
        std::array<T, PSO_DIM> point{};
        for (size_t i = 0; i < PSO_DIM; i++){ point[i] = x[i]; }
        return {point};
    }

    /// @brief One minimum is one row.
    size_t stride() const override{
        return 1;
    }

    /// @brief Maps a name onto a function; throws when it is not one of them.
    static Function parse(const std::string& name);

private:
    Function which;
};

extern template class AnalyticObjective<Real>;
extern template class AnalyticObjective<LongReal>;

#endif // PSO_OBJECTIVE_ANALYTIC_HPP
