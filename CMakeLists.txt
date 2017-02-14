cmake_minimum_required(VERSION 3.6)
project(lbfgsb_cpp_wrapper CXX Fortran)

enable_language(Fortran)
set(CMAKE_CXX_STANDARD 11)

# Set armadillo
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(SOURCE_FILES l_bfgs_b_example.cpp include/problem.h Lbfgsb.3.0/blas.f Lbfgsb.3.0/lbfgsb.f Lbfgsb.3.0/linpack.f Lbfgsb.3.0/timer.f
                 l_bfgs_b.h tests/src/test_functions.h)

add_executable(example ${SOURCE_FILES})
target_link_libraries(example ${ARMADILLO_LIBRARIES})
