cmake_minimum_required(VERSION 3.6)

set(SIMPLE_SOURCE_FILES simple_example.cpp)
add_executable(simple_example ${SIMPLE_SOURCE_FILES})
target_link_libraries(simple_example ${PROJECT_NAME})