cmake_minimum_required(VERSION 3.6)
project(lbfgsb_cpp_wrapper_example CXX Fortran)

enable_language(Fortran)
set(CMAKE_CXX_STANDARD 11)

file(GLOB FORTRAN_SRC
        "../Lbfgsb.3.0/*.f"
        )
set(SIMPLE_SOURCE_FILES simple_example.cpp ${FORTRAN_SRC})
add_executable(simple_example ${SIMPLE_SOURCE_FILES})
