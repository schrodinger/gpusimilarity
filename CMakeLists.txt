cmake_minimum_required(VERSION 3.11.1 FATAL_ERROR)
project(FastSimProject VERSION 1.0.0.0 DESCRIPTION "fingerpring similarity searching" LANGUAGES CXX CUDA)
include(CTest)
include(GNUInstallDirs)

set(QT_MINIMAL_VERSION 5.9.5)
set(BOOST_MINIMAL_VERSION 1.62.0)
set(CUDA_MINIMAL_VERSION 8.0.61)


# ccmake -DCMAKE_CUDA_COMPILER=path/to/nvcc ../gpusimilarity
# to also build unit test add:
# -DBOOST_ROOT=path/to/boost
# to generate the doc add:
# -DDOXYGEN_EXECUTABLE=/path/to/doxygen
# make doc_doxygen

if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS ${CUDA_MINIMAL_VERSION})
  message(FATAL_ERROR "Minimum supported CUDA version is ${CUDA_MINIMAL_VERSION}!")
endif()

# Qt libraries
find_package(Qt5Network REQUIRED)
if (Qt5Network_FOUND)
    if (Qt5Network_VERSION VERSION_LESS ${QT_MINIMAL_VERSION})
        message(FATAL_ERROR "Minimum supported Qt5 version is ${QT_MINIMAL_VERSION}!")
    endif()
else()
    message(SEND_ERROR "The Qt5Network library could not be found!")
endif(Qt5Network_FOUND)

# Workaround issue in Qt which affects Cuda
# ref: https://gitlab.kitware.com/cmake/cmake/issues/16915
# convert the fpic option in Qt5::Core over to
# INTERFACE_POSITION_INDEPENDENT_CODE
get_property(core_options TARGET Qt5::Core PROPERTY INTERFACE_COMPILE_OPTIONS)
string(REPLACE "-fPIC" "" new_core_options ${core_options})
set_property(TARGET Qt5::Core PROPERTY INTERFACE_COMPILE_OPTIONS ${new_core_options})
set_property(TARGET Qt5::Core PROPERTY INTERFACE_POSITION_INDEPENDENT_CODE "ON")

find_package(Doxygen)
if (DOXYGEN_FOUND) # make doc_doxygen
  message("Doxygen found!")
  # .in file generated with doxygen -g Doxyfile.in and later
  # modified to only include parameters with non-default values
  set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
  set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
  add_custom_target(doc_doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
endif()

add_library(fastsim fastsim.cpp fingerprintdb_cuda.cu)
target_include_directories(fastsim PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>)


add_executable(fastsim_exe main.cpp)
# Force C++11
target_compile_features(fastsim_exe PRIVATE cxx_override)

target_link_libraries(fastsim Qt5::Network)
target_link_libraries(fastsim_exe fastsim)

install(TARGETS fastsim_exe RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS fastsim EXPORT fastsimConfig
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES fastsim.h types.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT fastsimConfig DESTINATION share/fastsim/cmake)
export(TARGETS fastsim FILE fastsimConfig.cmake)

option(BUILD_ON_CHROOT_CENTOS6 "" false)
if(BUILD_ON_CHROOT_CENTOS6)
  install(FILES /usr/lib64/libicui18n.so.42 /usr/lib64/libicui18n.so.42.1 /usr/lib64/libicuuc.so.42 /usr/lib64/libicuuc.so.42.1 /usr/lib64/libicudata.so.42 /usr/lib64/libicudata.so.42.1 DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
