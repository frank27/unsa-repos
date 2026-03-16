cmake_minimum_required(VERSION 3.5)

project(MultiplyLapack LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find LAPACK (automatically brings in BLAS as a dependency)

find_package(LAPACK REQUIRED)

# Compiler flags  (-O2 -Wall equivalent for all build types)
add_compile_options(-O2 -Wall)

set(project_sources
    src/main.cpp
    src/DGEMMRoutine.h
    src/StandardAlgorithm.h
    src/Utils.h
)

add_executable(MultiplyLapack ${project_sources})

# Link LAPACK (and its BLAS dependency) to the executable
target_link_libraries(MultiplyLapack PRIVATE LAPACK::LAPACK)

target_compile_definitions(MultiplyLapack PRIVATE
    DATA_DIR="${CMAKE_SOURCE_DIR}/src/"
)