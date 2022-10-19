# Main CMakeLists file for OpenCAEPoro
#
# Sample usages:
#   cmake .                             // build in default configuration
#   cmake -DCMAKE_BUILD_TYPE=Release .  // build in Release configuration
#   cmake -DCMAKE_BUILD_TYPE=Debug .    // build in Debug configuration
#   cmake -DCMAKE_VERBOSE_MAKEFILE=ON . // build with verbose on
#   cmake -DUSE_FASP4BLKOIL=ON .        // build with FASP4BLKOIL support
#   cmake -DUSE_FASP4CUDA=ON .          // build with FASP4CUDA support
#   cmake -DUSE_UMFPACK=ON .            // build with UMFPACK support
#   cmake -DUSE_OPENMP=ON .             // build with OpenMP support

###############################################################################
## General environment setting
###############################################################################

# Minimum cmake version needed
cmake_minimum_required(VERSION 3.13)
cmake_policy(SET CMP0076 NEW)

# Name of the project: should appear after find_program
project(OpenCAEPoro C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)

# Set module lookup path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/modules)

###############################################################################
## External dependencies
###############################################################################

# Find BLAS and LAPACK implementations
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Add math and stdc++ libraries when not included by default
if (UNIX OR LINUX)
	set(ADD_STDLIBS m stdc++)
endif()

# Install external packages (faspsolver and ecl)
add_subdirectory(external)

# Other optional external packages 
include(optional)

###############################################################################
## Find compiler support for OpenMP. Since cmake 3.9, we can use 
## find_package(OpenMP). Keep this for older cmake!!!
###############################################################################

if(USE_OPENMP)
    if(CMAKE_VERSION VERSION_GREATER 3)
        cmake_policy(SET "CMP0054" NEW)
    endif()

    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")	 
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    else(OPENMP_FOUND)
        message("   WARNING: OpenMP was requested but disabled!")
    endif(OPENMP_FOUND)
endif(USE_OPENMP)

###############################################################################
## Set configuration types and default compiler flags
###############################################################################

set(CMAKE_CONFIGURATION_TYPES "Release;Debug" CACHE STRING "Configs" FORCE)
if (MSVC)
    # Compiler flags for Windows MSVC
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O3")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /DDEBUG /W4")
else()
    # Compiler flags for Linux or Mac OS X
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -Wall -g")
endif()

###############################################################################
## Project specific parameters
###############################################################################

# OpenCAEPoro library targets
add_library(OpenCAEPoro STATIC)
if(USE_FASP4CUDA)
    target_link_libraries(OpenCAEPoro PUBLIC fasp cublas cusparse cudart cudadevrt)
else()
    target_link_libraries(OpenCAEPoro PUBLIC fasp)
endif()
if(BUILD_ECL)
    target_link_libraries(OpenCAEPoro PUBLIC ecl)
endif()
target_link_libraries(OpenCAEPoro PUBLIC ${VTK_LIBRARIES})

add_subdirectory(src)
add_subdirectory(include)
add_subdirectory(main)

# set default build type: Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "the build type" FORCE)
endif()
message(STATUS "Set build type to ${CMAKE_BUILD_TYPE} configuration")

option(BUILD_TEST "Build test" OFF)
if(BUILD_TEST)
    include(CTest)
    # add_subdirectory(test)
endif()

###############################################################################
## Installtion targets for lib and executable files
###############################################################################

# Set install location for OpenCAEPoro library
set(CMAKE_INSTALL_LIBDIR "${PROJECT_SOURCE_DIR}/lib/")
install(TARGETS OpenCAEPoro
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# CMake uninstall target
if(NOT TARGET uninstall)
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/modules/cmake_uninstall.cmake.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
                   IMMEDIATE @ONLY)
    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
