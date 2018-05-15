cmake_minimum_required(VERSION 3.11.1 FATAL_ERROR)
project(FastSimTestProject VERSION 1.0.0.0 DESCRIPTION "unit test" LANGUAGES CXX)

find_package(Boost ${BOOST_MINIMAL_VERSION} COMPONENTS unit_test_framework REQUIRED)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")
include(BoostUnitTest)
include_directories(${Boost_INCLUDE_DIRS})
add_boost_test(test_fastsim.cpp fastsim ${Boost_LIBRARY_DIRS}:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

configure_file(small.fsim small.fsim COPYONLY)
