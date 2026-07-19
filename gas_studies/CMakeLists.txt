cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(p2_gas_scan CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Provided by the LCG view sourced in run_lxplus.sh (Garfield::Garfield target).
find_package(Garfield REQUIRED)

add_executable(p2_gas_scan p2_gas_scan.cpp)
target_link_libraries(p2_gas_scan Garfield::Garfield)
