cmake_minimum_required(VERSION 3.20)
project(PSO LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Definitions
add_compile_definitions(MPREAL NO_ODE_WARN)

# Platform-specific settings
if(APPLE)
    add_compile_definitions(__MAC__)
    include_directories(/opt/homebrew/Cellar/ /opt/homebrew/Cellar/mpfr/4.2.2/include/)
    link_directories(/usr/local/lib)
else()
    find_package(OpenMP REQUIRED)
endif()

set(ODEPACK_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/external/odepack)
set(LAZY_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/external/xdiff/external/lazy)
set(EIGEN_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen)

set(ODEPACK_ENABLE_MPREAL ON CACHE BOOL "Enable MPFR multi-precision support" FORCE)
add_subdirectory(${ODEPACK_ROOT} external/odepack)

add_library(lazy INTERFACE)
add_library(lazy::lazy ALIAS lazy)
target_include_directories(lazy INTERFACE ${LAZY_ROOT}/include)
target_compile_definitions(lazy INTERFACE LAZY_MPFR_RND=MPFR_RNDN)

add_library(eigen INTERFACE)
add_library(eigen::eigen ALIAS eigen)
target_include_directories(eigen INTERFACE ${EIGEN_ROOT})

# Release/Debug flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -funroll-loops -fno-math-errno")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fno-omit-frame-pointer")

# Source files
set(FUNCTIONS_SRC src/functions/functions.cpp)
set(GRID_SEARCH_SRC src/grid_search/grid_search.cpp)
set(GLOBAL_CLASSIC_SRC src/pso/global_classic/global_classic.cpp)
set(LOCAL_CLASSIC_SRC src/pso/enhanced/local_classic/local_classic.cpp)
set(DEFLECTION_SRC src/pso/enhanced/deflection/deflection.cpp)

# Libraries
set(MPFR_LIBS mpfr gmp)

add_library(pso_deps INTERFACE)
add_library(pso::deps ALIAS pso_deps)
target_link_libraries(pso_deps INTERFACE odepack::odepack lazy::lazy eigen::eigen ${MPFR_LIBS})
target_include_directories(pso_deps INTERFACE ${ODEPACK_ROOT}/include)

# def_demo
add_executable(def_demo
    demos/deflection_demo.cpp
    ${FUNCTIONS_SRC}
    ${DEFLECTION_SRC}
)
target_link_libraries(def_demo pso::deps)
if(NOT APPLE)
    target_link_libraries(def_demo OpenMP::OpenMP_CXX)
endif()

# grid_demo
add_executable(grid_demo
    demos/grid_search_demo.cpp
    ${FUNCTIONS_SRC}
    ${GRID_SEARCH_SRC}
    ${GLOBAL_CLASSIC_SRC}
    ${LOCAL_CLASSIC_SRC}
)
target_link_libraries(grid_demo pso::deps)
if(NOT APPLE)
    target_link_libraries(grid_demo OpenMP::OpenMP_CXX)
endif()

# local_demo
add_executable(local_demo
    demos/local_classic_demo.cpp
    ${FUNCTIONS_SRC}
    ${LOCAL_CLASSIC_SRC}
)
target_link_libraries(local_demo pso::deps)
if(NOT APPLE)
    target_link_libraries(local_demo OpenMP::OpenMP_CXX)
endif()

# global_demo
add_executable(global_demo
    demos/global_classic_demo.cpp
    ${FUNCTIONS_SRC}
    ${GLOBAL_CLASSIC_SRC}
)
target_link_libraries(global_demo pso::deps)
if(NOT APPLE)
    target_link_libraries(global_demo OpenMP::OpenMP_CXX)
endif()
