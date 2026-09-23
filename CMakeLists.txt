cmake_minimum_required(VERSION 3.20)

project(pso LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Writes build/compile_commands.json, which clangd discovers on its own.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Matching the reference implementation's flags, so timings compare directly.
# -fno-math-errno lets sqrt inline to a single instruction, and the raised
# inline limit matters because the stepper is deep templates over small
# functions. Neither relaxes IEEE semantics, so results are unaffected.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -funroll-loops -fno-math-errno")
set(CMAKE_CXX_FLAGS_DEBUG   "-g -O0 -fno-omit-frame-pointer")

# The objective is evaluated one particle per thread, so the parallel-for in
# macOS: Apple's own clang ships no OpenMP runtime, and Homebrew installs into
# a prefix that differs by architecture (/opt/homebrew on Apple silicon,
# /usr/local on Intel) and by package version. Ask brew where things are
# rather than hardcoding a path that goes stale on the next upgrade.
if(APPLE)
    find_program(BREW_EXECUTABLE brew)

    if(BREW_EXECUTABLE)
        foreach(pkg libomp mpfr gmp)
            execute_process(
                COMMAND ${BREW_EXECUTABLE} --prefix ${pkg}
                OUTPUT_VARIABLE PSO_BREW_PREFIX
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
                RESULT_VARIABLE PSO_BREW_RESULT)

            if(PSO_BREW_RESULT EQUAL 0 AND EXISTS "${PSO_BREW_PREFIX}")
                list(APPEND CMAKE_PREFIX_PATH "${PSO_BREW_PREFIX}")
                include_directories(SYSTEM "${PSO_BREW_PREFIX}/include")
                link_directories("${PSO_BREW_PREFIX}/lib")
            endif()
        endforeach()
    endif()
endif()

# The objective is evaluated one particle per thread. OpenMP is the preferred
# way to do that, but it is genuinely optional on macOS, where Grand Central
# Dispatch is always available and needs no third-party runtime. Without
# either, the pragma is silently ignored and every run is serial -- so the
# fallback is explicit rather than accidental.
find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    set(PSO_PARALLEL_BACKEND "OpenMP")
elseif(APPLE)
    set(PSO_PARALLEL_BACKEND "Grand Central Dispatch")
else()
    message(FATAL_ERROR
        "OpenMP not found. It is required on this platform; on macOS the build "
        "falls back to Grand Central Dispatch, elsewhere there is no fallback. "
        "Install it (libomp-dev on Debian/Ubuntu) and reconfigure.")
endif()

message(STATUS "Objective evaluation parallelised with: ${PSO_PARALLEL_BACKEND}")

# The reference implementation defined NO_ODE_WARN; odecraft now spells that
# as an option. Warnings on every rejected step would otherwise dominate the
# console, since particles off the energy shell are routine here.
set(ODECRAFT_NO_WARN ON CACHE BOOL "" FORCE)

# Round-to-nearest for lazex's mpreal expressions, matching the original's
# LAZY_MPFR_RND. This is also lazex's own default; pinned so a change upstream
# cannot silently alter results.
set(LAZEX_MPFR_RND "MPFR_RNDN" CACHE STRING "" FORCE)

# The exact phase of the precision swap. ON gives LongReal =
# lazex::LazyType<mpfr::mpreal> at mpfr_prec bits and needs MPFR and GMP; OFF
# gives long double, which drops both dependencies but caps precision at 64
# bits of mantissa. Configure with -DPSO_USE_MPREAL_SWAP=OFF.

option(PSO_USE_MPREAL_SWAP "Use MPFR for the exact phase instead of long double" ON)

if(PSO_USE_MPREAL_SWAP)
    # mpreal.h reaches us through lazex::mpfr, and lazex is configured from
    # inside odecraft, so this has to be set before add_subdirectory below.
    set(LAZEX_ENABLE_MPREAL ON CACHE BOOL "" FORCE)
endif()

add_subdirectory(external/odecraft)
add_subdirectory(external/tomlplusplus)

# Everything that is not a template lives here, compiled once and shared by
# every executable. particle.hpp and swarm.hpp stay header-only: they are
# templates over the scalar type and the stepper, so there is nothing to
# compile separately without explicit instantiation.
add_library(pso_core STATIC
    src/input/parameters.cpp
    src/input/initialize.cpp
    src/input/validate.cpp
    src/objective/objective.cpp
    src/objective/poincare.cpp
    src/objective/analytic.cpp
    src/integrate/events.cpp
    src/integrate/make_integrator.cpp

    # One per stepper: each instantiates its stepper for both scalar types
    # with CalderaODE::Rhs inlined, and they compile in parallel.
    src/integrate/euler.cpp
    src/integrate/rk4.cpp
    src/integrate/rk23.cpp
    src/integrate/rk45.cpp
    src/integrate/dop853.cpp
    src/integrate/bdf.cpp

    src/swarm/swarm.cpp
    src/swarm/topology.cpp
    src/driver/search.cpp
    src/driver/grid.cpp
    src/driver/dispatch.cpp
)

# Both spellings resolve:
#   #include <src/input/parameters.hpp>   (repo root)
#   #include <input/parameters.hpp>       (src/)
target_include_directories(pso_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# PUBLIC throughout: consumers include these headers, and PSO_USE_MPREAL_SWAP
# must agree in every translation unit -- a mismatch would give LongReal two
# different meanings and quietly violate the ODR.
target_link_libraries(pso_core PUBLIC
    odecraft::odecraft
    tomlplusplus::tomlplusplus
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(pso_core PUBLIC OpenMP::OpenMP_CXX)
else()
    # Blocks and libdispatch are part of the platform; nothing to link.
    target_compile_definitions(pso_core PUBLIC PSO_USE_GCD)
endif()

target_compile_definitions(pso_core PUBLIC
    PSO_USE_MPREAL_SWAP=$<BOOL:${PSO_USE_MPREAL_SWAP}>
)

if(PSO_USE_MPREAL_SWAP)
    target_link_libraries(pso_core PUBLIC lazex::mpfr)
endif()

add_executable(pso
    src/main.cpp
)

target_link_libraries(pso PRIVATE pso_core)

# Analogue of the reference implementation's local_classic_demo, for timing.
# Demos never hardcode parameters: this one reads demos/local_demo.toml, whose
# path is baked in so the binary takes no arguments and works from any cwd.
# Demos are EXCLUDE_FROM_ALL: not part of the default build.
#   cmake --build build --target local_demo
# Each reads its own .toml, whose path is baked in so the binary takes no
# arguments and works from any working directory.
foreach(demo local global deflection grid)
    add_executable(${demo}_demo EXCLUDE_FROM_ALL demos/${demo}_demo.cpp)
    target_compile_definitions(${demo}_demo PRIVATE
        PSO_DEMO_INPUT="${CMAKE_CURRENT_SOURCE_DIR}/demos/${demo}_demo.toml"
    )
    target_link_libraries(${demo}_demo PRIVATE pso_core)
endforeach()
