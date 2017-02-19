cmake_minimum_required(VERSION 3.6)
project(lbfgsb_cpp_wrapper CXX Fortran)

enable_language(Fortran)
set(CMAKE_CXX_STANDARD 11)

# Set armadillo
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
# Set eigen
SET(EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
IF(NOT EIGEN3_INCLUDE_DIR)
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF()
INCLUDE_DIRECTORIES("${EIGEN3_INCLUDE_DIR}")

set(SOURCE_FILES l_bfgs_b_example.cpp include/problem.h Lbfgsb.3.0/blas.f Lbfgsb.3.0/lbfgsb.f Lbfgsb.3.0/linpack.f Lbfgsb.3.0/timer.f
        include/l_bfgs_b.h tests/src/test_functions.h)

add_executable(example ${SOURCE_FILES})
target_link_libraries(example ${ARMADILLO_LIBRARIES})
