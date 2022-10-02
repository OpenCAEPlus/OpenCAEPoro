# Main CMakeLists file for OpenCAEPoro
# 
# Author            Date           Action
# -----------------------------------------------------------------------------
# Shizhe Li         Sep/21/2021    Create CMakeLists
# Chensong Zhang    Sep/24/2021    Fix FASP linkage problems
# Chensong Zhang    Oct/16/2021    Fix duplications
# Chensong Zhang    Jan/19/2022    Use FASP open-source part as default solver
# Chensong Zhang    Feb/18/2022    Fix math.lib linkage problem
#
# Some sample usages:
#   cmake .                             // build in default configuration
#   cmake -DCMAKE_BUILD_TYPE=Release .  // build in Release configuration
#   cmake -DCMAKE_BUILD_TYPE=Debug .    // build in Debug configuration
#   cmake -DCMAKE_VERBOSE_MAKEFILE=ON . // build with verbose on
#   cmake -DUSE_FASP4BLKOIL=ON .        // build with FASP4BLKOIL support
#   cmake -DUSE_UMFPACK=ON .            // build with UMFPACK support
#   cmake -DUSE_OPENMP=ON .             // build with OpenMP support

###############################################################################
## General environment setting
###############################################################################

# Minimum cmake version needed
cmake_minimum_required(VERSION 3.13)
cmake_policy(SET CMP0076 NEW)
# Default C compiler look-up order
find_program(THE_C NAMES $ENV{CC} gcc clang icc icx icl)
set(CMAKE_C_COMPILER ${THE_C} CACHE INTERNAL "the C compiler" FORCE)

# Default C++ compiler look-up order
find_program(THE_CXX NAMES $ENV{CXX} g++ clang++ icpc icx icl)
set(CMAKE_CXX_COMPILER ${THE_CXX} CACHE INTERNAL "the CXX compiler" FORCE)

# Name of the project: should appear after find_program
project(OpenCAEPoro C CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)

# Set module lookup path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/modules)

###############################################################################
## External libraries
###############################################################################

# Find BLAS and LAPACK implementations
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)




###############################################################################
## Find compiler support for OpenMP. Since cmake 3.9, we can use 
## find_package(OpenMP REQUIRED). Keep this for older cmake!!!
###############################################################################

if(USE_OPENMP)
    if(CMAKE_VERSION VERSION_GREATER 3)
        cmake_policy(SET "CMP0054" NEW)
    endif()

    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")	 
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
    else(OPENMP_FOUND)
        message("   WARNING: OpenMP was requested but disabled!")
    endif(OPENMP_FOUND)
endif(USE_OPENMP)

###############################################################################
## Set configuration types and their default compiler flags
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

## requires math library
if (UNIX OR LINUX)
	set(ADD_STDLIBS m stdc++)
endif()

# Set default build type to Release
set(build_type "Release")
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "${build_type}" CACHE STRING "the build type" FORCE)
endif()
message(STATUS "Set build type to ${CMAKE_BUILD_TYPE} configuration")

###############################################################################
## Project specific parameters
###############################################################################

# The OpenCAEPoro library target
add_subdirectory(external)
add_library(OpenCAEPoro STATIC)
target_link_libraries(OpenCAEPoro PUBLIC fasp)
add_subdirectory(src)
add_subdirectory(include)
add_subdirectory(main)

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