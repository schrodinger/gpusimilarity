cmake_minimum_required(VERSION 3.10.2 FATAL_ERROR)
project(GPUSimTestProject VERSION 1.0.0.0 DESCRIPTION "unit test" LANGUAGES CXX)

find_package(Boost ${BOOST_MINIMAL_VERSION} COMPONENTS unit_test_framework REQUIRED)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}")
include(BoostUnitTest)
include_directories(${Boost_INCLUDE_DIRS})
add_boost_test(test_gpusim.cpp gpusim ${Boost_LIBRARY_DIRS}:${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

configure_file(small.fsim small.fsim COPYONLY)
# Force -fPIC on executable. This is Qt bug.
# https://bugreports.qt.io/browse/QTBUG-51593
# https://cmake.org/Bug/view.php?id=9659
# https://gitlab.kitware.com/cmake/cmake/issues/16002
if(Qt5_POSITION_INDEPENDENT_CODE)
  target_compile_options(test_gpusim PRIVATE -fPIC)
endif()
