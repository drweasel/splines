cmake_minimum_required(VERSION 3.10)
project(Splines LANGUAGES CXX VERSION 1.0.0 DESCRIPTION "Spline Approximation")

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(GTest)
enable_testing()

include_directories(${Splines_SOURCE_DIR}/include)

add_compile_options(-std=c++14)

if (${GTEST_FOUND})
  message("configuring unit tests ...")
  include_directories(${GTEST_INCLUDE_DIRS})
  add_subdirectory(test)
endif()

