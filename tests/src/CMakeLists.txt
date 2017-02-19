enable_testing()

file(GLOB FORTRAN_SRC
        "../../Lbfgsb.3.0/*.f"
)
set(SOURCE_TEST_FILES ${FORTRAN_SRC}
                      test_gtest.cpp)
add_executable(run_test ${SOURCE_TEST_FILES})
target_include_directories(run_test PUBLIC ${gtests_SOURCE_DIR})
target_include_directories(run_test PUBLIC ${ARMADILLO_INCLUDE_DIRS})
target_include_directories(run_test PUBLIC ${EIGEN3_INCLUDE_DIR})
target_link_libraries(run_test gtest_main ${ARMADILLO_LIBRARIES})

