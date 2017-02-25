cmake_minimum_required(VERSION 3.6)

IF (BUILD_SIMPLE_EX)
    add_executable(simple_example simple_example.cpp)
    target_link_libraries(simple_example ${PROJECT_NAME})
ENDIF()

IF (BUILD_FULL_EX)
    add_executable(full_example full_example.cpp)
    target_include_directories(full_example PUBLIC ${ARMADILLO_INCLUDE_DIRS})
    target_include_directories(full_example PUBLIC ${EIGEN3_INCLUDE_DIR})
    target_link_libraries(full_example ${PROJECT_NAME})
    target_link_libraries(full_example ${ARMADILLO_LIBRARIES})
ENDIF()

