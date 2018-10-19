cmake_minimum_required(VERSION 3.10)

include_directories(${CMAKE_SOURCE_DIR})

add_executable(unit_tests
  natural_spline_test.cpp
  bezier_spline_test.cpp
  polynomials_test.cpp
  zeros_test.cpp
  )

target_link_libraries(unit_tests ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES})

