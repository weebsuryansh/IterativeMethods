cmake_minimum_required(VERSION 3.29)
project(IterativeMethods)

set(CMAKE_CXX_STANDARD 20)

include_directories(${CMAKE_SOURCE_DIR}/Eigen)

add_executable(IterativeMethods
        main.cpp
        jacobi.cpp
        gauss_seidel.cpp
        conjugate_gradient.cpp
        gmres.cpp
)
