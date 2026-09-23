# Particle Swarm Optimisation [![status: active](https://github.com/GIScience/badges/raw/master/status/active.svg)](https://github.com/GIScience/badges#active) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Particle swarm optimisation, applied to finding **periodic orbits** in a
two-degree-of-freedom Hamiltonian system.

## Overview

A swarm of candidate points wanders a 2D search space. Each one is scored by a
function you choose, the good ones attract the others, and the swarm converges
on a minimum. That is particle swarm optimisation, and it is about fifty lines
of arithmetic.

Everything else in this repository exists because of *which* function is being
minimised.

The default one works like this. Pick a starting point in a physical system.
Let it evolve in time. Wait until it comes back to a chosen slice of the space
— a *surface of section* — and measure how far it landed from where it began.
If that distance is zero, the trajectory closed on itself: you have found a
**periodic orbit**. So "minimise the return distance" and "find a periodic
orbit" are the same instruction.

Two consequences shape the whole program:

- **Evaluating the function is expensive.** Each one integrates a differential
  equation. A run does hundreds of thousands of them, so the objective, not the
  swarm, is where all the time goes.
- **The answer must be very precise.** `double` runs out of accuracy around
  1e-6, and the target is 1e-8 or better. The program therefore starts in
  `double` and switches to arbitrary-precision arithmetic partway through.
  That switch is the most unusual thing here and is explained below.

## Quick start

You need a C++20 compiler, CMake 3.20 or newer, and the MPFR and GMP
development libraries. For Debian / Ubuntu:

```sh
sudo apt install libmpfr-dev libgmp-dev libomp-dev      # Debian / Ubuntu
```

For macOS:
```sh
brew install mpfr gmp
```

Build the project:
```sh
git clone --recursive https://github.com/phyzan/Particle_Swarm_Optimization
cd Particle_Swarm_Optimization
cmake -S . -B build
cmake --build build -j4
```
and run the executable from the example parameter file `input.toml`:
```sh
./build/pso input.toml
```
or build and run a demo:
```sh
cmake --build build --target local_demo
./build/local_demo
```


Configure with much greater performance using
```sh
cmake -S . -B build -DPSO_USE_MPREAL_SWAP=OFF
```
as explained [later](#The-precision-swap).


On macOS, `libomp` is optional: Apple's clang ships no OpenMP runtime, so the
build falls back to Grand Central Dispatch, which is always present. Install
`libomp` as well if you would rather use OpenMP. Configure prints which one it
picked:

```
-- Objective evaluation parallelised with: Grand Central Dispatch
```

On that path `execution.threads` is advisory rather than a cap — libdispatch
sizes its own pool from the machine and the current load.

### Reading the output

```
  [double] residual 4.83e-08 below swap_tol after 202 iterations; continuing at 64 bits
  [mpreal] converged after 214 iterations, residual 8.51e-09

1 minimum found
  orbit 1
    x = -0.000128774932   px = -0.000102234032
```

- **residual** — how far the current best candidate is from solving the
  problem. Smaller is better; the run stops when it drops under `objfun_tol`.
- **`[double]` then `[mpreal]`** — the two precision phases. The handover
  happened at iteration 202; the iteration count carries across, so the run
  used 214 in total.
- **orbit** — the answer. For this problem a point is `(x, px)`: a position and
  a momentum.

## The input file

Everything is set in a TOML file. Nothing is hardcoded, and the program reads
exactly one file, named on the command line (default `input.toml`).

```toml
[domain]                      # where to search
lower = [-2.5, -0.5]          # required: no sensible default exists
upper = [ 2.5,  0.5]

[convergence]                 # when to stop
objfun_target = 0.0           # the value to reach
objfun_tol    = 1e-8          # success when |value - target| <= this
swap_tol      = 1e-5          # switch to high precision below this
mpfr_prec     = 100           # bits of precision after the switch
max_iter      = 15000         # budget, shared by BOTH phases

[dynamics]                    # how the swarm moves
pop_size = 20                 # number of candidate points
c1 = 2.0                      # pull toward a point's own best find
c2 = 1.7                      # pull toward its neighbours' best find

[strategy]
kind = "single"               # single | deflection | grid

[objective]
kind = "poincare"             # poincare | analytic

[poincare]                    # the physics
energy      = 17.0
section     = -1.8019693
n_crossings = 1               # 1 = simple orbit, 2 = orbit of double the period
integrator  = "RK45"          # Euler RK4 RK23 RK45 DOP853 BDF
rtol        = 1e-12           # integration accuracy

[execution]
seed    = 42                  # runs are reproducible
threads = 8
```

Every field has a default except `domain.lower` and `domain.upper`. The
defaults, their allowed ranges and what each one does are documented next to
the field in [`src/input/parameters.hpp`](src/input/parameters.hpp), which is
the authoritative list.

**A misspelled name is an error, not a silent default:**

```
error: unrecognised parameters in input.toml:
  dynamics.c11   (line 14)
  dynmaics.c2    (line 22)
```

Contradictory settings are caught before anything runs:

```
error: convergence.swap_tol must exceed objfun_tol, else the run converges
       before it can swap
```

## Choosing what to minimise

```toml
[objective]
kind = "poincare"
```

**`poincare`** — the real problem. Integrates the trajectory and measures the
return distance. Expensive, and what the physics is about.

**`analytic`** — three standard test functions with known answers, evaluated in
a few arithmetic operations instead of a whole integration:

```toml
[objective]
kind = "analytic"
[analytic]
function = "rastrigin"        # sphere | rastrigin | rosenbrock
```

Use these to watch the swarm work without waiting on the physics. A run
finishes in milliseconds, so they are the fastest way to get a feel for what
the parameters do. All three have their minimum at value 0, so
`objfun_target = 0` suits them all.

## Choosing how to search

```toml
[strategy]
kind = "single"
```

**`single`** — one swarm, one minimum. Start here.

**`deflection`** — runs the swarm repeatedly. After each success it makes the
function blow up near the minimum it already found, so the next run cannot
converge to the same place. Finds several minima in one go.

```toml
[strategy]
kind = "deflection"
[deflection]
runs = 3
```

**`grid`** — chops the domain into a grid and runs an independent swarm in each
cell. The most thorough option and by far the slowest.

```toml
[strategy]
kind = "grid"
[grid]
resolution = 2                # 2^(2*resolution) cells: 2 -> 16, 4 -> 256
batch      = 0                # cells at once; 0 = fill the machine
```

`resolution` is an exponent, so it grows fast: each step up is four times the
work in 2D. Start at 2.

Each cell writes its own log to `execution.output_dir`, so parallel runs never
interleave.

## The precision swap

This is the idea worth understanding.

`double` is fast but only carries about 16 digits, and its own rounding error
swamps the answer below roughly 1e-6. Arbitrary-precision arithmetic reaches
any accuracy you ask for but costs perhaps 200 times more per operation.

So the program uses both. It runs in `double` until the residual falls below
`swap_tol`, hands the **entire swarm** — positions, velocities, every point's
best-known result — to a high-precision copy, and carries on. The iteration
count crosses with it, which is why `max_iter` is one shared budget rather than
two.

```
double, fast          swap_tol           mpreal, exact
├──────────────────────────┤──────────────────────────┤
start                    1e-5                   objfun_tol
```

Two settings worth knowing:

- `swap_tol` at or below `objfun_target` never swaps — a pure `double` run.
- `swap_tol` above the starting residual swaps immediately — everything slow
  and exact.

If you do not need more than about 19 digits, you can drop the MPFR dependency
entirely and use the hardware's widest float:

```sh
cmake -S . -B build -DPSO_USE_MPREAL_SWAP=OFF
```

That build is roughly ten times faster and links neither MPFR nor GMP. It caps
precision at 64 bits of mantissa, which is ample for `objfun_tol = 1e-8`.

## Demos

Four ready-made configurations, each mirroring one of the original
implementation's examples. They are **not** built by default, because one of
them is a 256-cell grid run:

```sh
cmake --build build --target local_demo
./build/local_demo
```

| target | what it shows |
|---|---|
| `local_demo` | neighbourhood topology (LSH) |
| `global_demo` | one attractor for the whole swarm |
| `deflection_demo` | repeated runs, several minima |
| `grid_demo` | one swarm per subspace — **slow**, 256 cells |

Each takes no arguments and reads the `.toml` beside it, whose path is compiled
in. The `.toml` files document which setting came from which line of the
original code, so they double as a migration record.

## Layout

| path | holds |
|---|---|
| `src/input/` | every parameter and its default; the TOML reader; validation |
| `src/objective/` | what is being minimised — the interface, Poincaré, analytic |
| `src/integrate/` | the ODE system, the section event, one file per stepper |
| `src/swarm/` | swarm state, the iteration loop, neighbourhood topologies |
| `src/driver/` | the precision swap; the single, deflection and grid strategies |
| `demos/` | four example configurations |

The layering is strict: the swarm and the drivers know nothing about physics.
They ask an objective for a value at a point, and that is the entire contract.

## Adding your own objective

The swarm needs three things from an objective, and nothing else:

```cpp
template<typename T>
class Objective{
    virtual T evaluate(const T* x) = 0;                 // value at x, +inf if infeasible
    virtual std::vector<std::array<T, PSO_DIM>> refine(const T* x) = 0;
    virtual size_t stride() const = 0;                  // rows per recorded minimum
};
```

To add one:

1. Write a class implementing those three, in `src/objective/`.
2. Add a value to `ObjectiveKind` in `src/input/parameters.hpp`, with a struct
   for its own parameters.
3. Read that struct in `read_parameters` and check it in `validate`.
4. Add one line to `make_objective` in `src/objective/objective.cpp`.

Nothing in `src/swarm/`, `src/driver/` or `src/integrate/` changes.
[`src/objective/analytic.hpp`](src/objective/analytic.hpp) is the worked
example and is short enough to read in a minute.

## Exit codes

| code | meaning |
|---|---|
| 0 | at least one minimum reached `objfun_tol` |
| 1 | budget exhausted or stalled; nothing converged |
| 2 | invalid input; nothing ran |

## Credits

Method and test problem from Katsanikas, Bakos and Wiggins (2026), *The
Computation of Periodic Orbits in Hamiltonian Systems Using Swarm
Intelligence*, International Journal of Bifurcation and Chaos 36, 2650102.

Integration by [odecraft](https://github.com/phyzan/odecraft). Input parsing by
[toml++](https://github.com/marzer/tomlplusplus).
