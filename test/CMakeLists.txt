enable_testing()
add_executable(
    Utils_test
    Utils_test.cpp
)
target_include_directories(Utils_test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)

target_link_libraries(Utils_test PRIVATE Eigen3::Eigen GTest::gtest_main RadiationLib)
include(GoogleTest)
gtest_discover_tests(Utils_test)
