# Particle Swarm Optimisation [![status: active](https://github.com/GIScience/badges/raw/master/status/active.svg)](https://github.com/GIScience/badges#active) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This is an implementation of the Particle Swarm Optimisation (PSO) algorithm, which was used to conduct the executions for the publication [publication].

The PSO algorithm is a Swarm intelligence method that solves global optimisation problems. More information about how the algorithm functions can be found in the [publication].

## Objective function

The algorithm converges to a minumum by utilizing the objective function provided by the user. The objective function is used to improve the position of the particles in the space.

As a demo objective function, a minimization function in a Poincaré Surface Section of a 4-dimensional phase-space is provided, the same as the one that is used in the calculations of the publication executions. The user may add their own objective function. There is guidance provided in the comments in `src/objective_functions` template file.

It is important to note that the parameters of PSO and the Objective function have to be correctly tuned in order for the algorithm to correctly converge and identify global (or local) minima. If the parameters are not set correctly the algorithm may have undefined behaviours like not converging to a minimum or converging only to one minimum, unable to further explore the defined space.

PSO may also work with objective functions that locate maxima instead of minima in a space.

## Algorithm variants

A small explaination of the PSO variants that are provided in this implementation.

### Classic Global PSO

The global variant of PSO is the classic variant of the algorithm as described in [Recent approaches to global optimization problems
through Particle Swarm Optimization](https://www.researchgate.net/publication/228746170_Recent_approaches_to_global_optimization_problems_through_Particle_Swarm_Optimization).

### Classic Local PSO

The local variant of the algorithm creates smaller neighbourhoods, resulting in more efficient searching of the space. This way if one neighbourhood gets stuck, the others can compensate by continuing the exploration of the space. To create said neighbourhoods, the LSH algorithm is used. Local PSO is also more optimal than the Global version in high-dimensional spaces, since it projects them to lower dimensions with the use of LSH, significantly improving performance.

### Deflection (with or without repulsion)

Deflection (and repulsion) is a technique that enhances the performance of Global or Local PSO and allows it to locate more than one minima. Essentially it is a measure that prevents the algorithm from converging to the same minimum multiple times. This technique is not as good as Grid Search for multiple minima discovery in a space.

### Grid Search

A method that disects the space into smaller subspaces and executes the global or local PSO variants in the latter. It is a more consistent way of finding minima in the space, if the space is disected into small enough subspaces. The way the number of subspaces is calculated is  $subspaces = 2^{dim * res}$, where `dim` is the dimension of the space and `res` is the resolution, a user defined variable. For example, if the space is 3D and the user defines resolution as 2, we have $2^{3*2}=64$ subspaces.

## Parameters

The parameters can be defined using the structure templates found in the header files. Most parameters have default values that might or might not be good for other objective functions.

- PSO parameters

| Parameter |                          Information                          |
|-----------|---------------------------------------------------------------|
|  dim      | The number of dimensions that the PSO particles have          |
|  max_it   | The maximum number of iterations of the algorithm             |
|  popsize  | The size of the population (Number of particles)              |
|  c1       | [Recent approaches to global optimization problems through Particle Swarm Optimization](https://www.researchgate.net/publication/228746170_Recent_approaches_to_global_optimization_problems_through_Particle_Swarm_Optimization)                                                 |
|  c2       | [Recent approaches to global optimization problems through Particle Swarm Optimization](https://www.researchgate.net/publication/228746170_Recent_approaches_to_global_optimization_problems_through_Particle_Swarm_Optimization)                                                 |
|  gm       | The minimum (or maximum) value that has to be achieved        |
|  max_w    | [Recent approaches to global optimization problems through Particle Swarm Optimization](https://www.researchgate.net/publication/228746170_Recent_approaches_to_global_optimization_problems_through_Particle_Swarm_Optimization)                                                 |
|  min_w    | [Recent approaches to global optimization problems through Particle Swarm Optimization](https://www.researchgate.net/publication/228746170_Recent_approaches_to_global_optimization_problems_through_Particle_Swarm_Optimization)                                                 |
|  err_goal | The error tolerance that determines if convergence is achieved|
|  bounds   | The boundaries of the space of the particles                  |

- Enhanced Parameters

| Parameter |                           Information                            |
|-----------|------------------------------------------------------------------|
|lsh_k      | The number of h() functions used for the hash of a hash-table    |
|lsh_w      | [Locality-Sensitive Hashing for Finding Nearest Neighbors](https://www.slaney.org/malcolm/yahoo/Slaney2008-LSHTutorial.pdf)                                                    |
|lsh_L      | The number of hash-tables                                        |
|con_k      | The constriction value for the velocity of PSO                   |
|rep_rho    | How strong the repulsion from the repulsion zone will be         |
|rep_radius | The hyper-sphere radius of the repulsion zone                    |
|lsh_radius | The hyper-sphere radius of the LSH search-space around a particle|

- Grid Parameters

| Parameter       |                                    Information                                     |
|-----------------|------------------------------------------------------------------------------------|
|resolution       | The resolution defines how much the space will be disected                         |
|batch_size       | How many threads will run at the same time                                         |
|starting_subspace| Sets the subspace that gridsearch will visit first, ignoring all the previous ones |
|dest             | The destination folder for the output of the program                               |
|type             | The algorithm variant that will be used, either `Global_Classic` or `Local_Classic`|
|gp               | A struct with PSO, Enhanced and Objective Function parameters                      |

- Objective Function Parameters
    - Poincare Parameters

| Parameter |                           Information                           |
|-----------|-----------------------------------------------------------------|
|  p        | Number of sections of an orbit and the PSS                      |
|  threads  | Number of threads used by the objective function (only on Linux)|
|  c1       |Controls parabolic curvature                                     |
|  c2       |Controls a linear tilt or slope in the y-direction               |
|  c3       |Controls higher-order distortion                                 |
|  dt       |The initial stepsize for the ODE integration                     |
|  ene      | Energy of a 2D caldera like potential energy surface (PES)      |
|  xpoin    | The Y coordinate                                                |
|  err_goal | The error tolerance that determines if convergence is achieved  |

For a more deep understanding of what each parameter does, it is advised to consult the [publication]. The default values are not suggested values, they are placeholders. Every system has different values that are efficient and they could differ a lot from the default ones.

- Swap Point

All methods have a user-set parameter called `swap_point`. Swap point is a variable swap from `double` to `MPReal`. `Double` variables are fast when used in calculations, but can accumulate a significant amount of error when the calculations require great accuracy ($1e-6$ or smaller).

On the other hand, `MPReal` variables make the calculations really slow (depending on the `precision` set by the user, generally around 200 times slower), but the calculations become 100% accurate.

Since both accuracy and speed are important, `swap_point` is essentially the point where `double` accuracy starts degrading. A good `swap_point` is $1e-6$. As a result, we can utilize both the speed of `double` variables and the accuracy of `MPReal` variables, without significantly sacrificing either of the two.

If the user wants to use only `double` variables, they can set the `swap_point` to be equal to the convergence point (e.g. if the algorithm converges to -1, then `swap_point = -1`). If the user wants to use only `MPReal` variables, set the swap_point to a significantly large number (e.g. 1000).

# Compiling the project

To compile the demo version, execute `make all` in the `/bin` folder if the environment is Linux. If it is MACOS, execute `make __MAC__=true all`. The difference between these two compilations is the implementation of the threads. In a Linux environment the user can define the number of threads they want to use, but in a MACOS environment the system itself regularises the number of threads used by the program.

# Testing and Production Environment

The testing and production of this project was done on Linux Ubuntu 22.04 and on MACOS version 15.5. The program was compiled using g++ 11.4.0 in Linux and the native clang++ compiler in MACOS (version 17.0.0 as of the publication of the project).

# Credits

- [Publication]
- [Foivos Zanias](https://github.com/phyzan), for the creation of the [ODE pack](https://github.com/phyzan/odepack) used for the implementation of this project, as well as insights for the general implementation.

# Versions

### Version 0.0.0
- Initial development version

### Version 0.1.0
- Added an option to start grid search from a selected subspace, ignoring the previous ones.

### Version 0.1.1
- Refactored code to follow Microsoft coding standards.

### Version 0.1.2
- Updated ReadMe and did minor code changes in grid_search.