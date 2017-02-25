cmake_minimum_required(VERSION 3.6)
project(lbfgsb_cpp CXX Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

option(BUILD_TESTS "Build tests" OFF)
option(BUILD_SIMPLE_EX "Build the simple example" ON)
option(BUILD_FULL_EX "Build the full example" OFF)

enable_language(Fortran)
set(CMAKE_CXX_STANDARD 11)

IF(BUILD_TESTS OR BUILD_FULL_EX)
    # Set armadillo
    find_package(Armadillo REQUIRED)
    # Set eigen
    SET(EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
    IF(NOT EIGEN3_INCLUDE_DIR)
        MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
    ENDIF()
ENDIF()


file(GLOB FORTRAN_SRC "Lbfgsb.3.0/*.f" )
file(GLOB HEADERS include/${PROJECT_NAME}/*.h)
set(SOURCE_FILES ${FORTRAN_SRC})

include_directories(include)
add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES})

# Install library
install(TARGETS ${PROJECT_NAME} DESTINATION lib/${PROJECT_NAME})
# Install library headers
install(FILES ${HEADERS} DESTINATION include/${PROJECT_NAME})

IF (BUILD_TESTS)
    add_subdirectory(tests)
ENDIF()

IF (BUILD_SIMPLE_EX OR BUILD_FULL_EX)
    add_subdirectory(examples)
ENDIF()