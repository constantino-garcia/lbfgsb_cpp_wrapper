cmake_minimum_required(VERSION 3.6)
project(lbfgsb_cpp_wrapper CXX Fortran)

enable_language(Fortran)
set(CMAKE_CXX_STANDARD 11)

# Set armadillo
find_package(Armadillo REQUIRED)
# Set eigen
SET(EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
IF(NOT EIGEN3_INCLUDE_DIR)
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF()

file(GLOB FORTRAN_SRC
        "Lbfgsb.3.0/*.f"
        )
set(FULL_SOURCE_FILES full_example.cpp ${FORTRAN_SRC})

add_executable(full_example ${FULL_SOURCE_FILES})
target_include_directories(full_example PUBLIC ${ARMADILLO_INCLUDE_DIRS})
target_include_directories(full_example PUBLIC ${EIGEN3_INCLUDE_DIR})
target_link_libraries(full_example ${ARMADILLO_LIBRARIES})

# add tests
add_subdirectory(tests)
