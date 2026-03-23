cmake_minimum_required(VERSION 3.5)

project(MultiplyAnalysis LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(LAPACK REQUIRED)

add_compile_options(-O2 -Wall)

set (project_sources
    src/main.cpp
    src/DiagProver.h
    src/DiagProver.cpp
    src/Diagonalization.h
    src/Diagonalization.cpp
)

add_executable(MultiplyAnalysis ${project_sources})

# Link LAPACK (and its BLAS dependency) to the executable
target_link_libraries(MultiplyAnalysis PRIVATE LAPACK::LAPACK)

target_compile_definitions(MultiplyAnalysis PRIVATE
    DATA_DIR="${CMAKE_SOURCE_DIR}/src/"
)