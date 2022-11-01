# ##############################################################################
# Find compiler support for OpenMP. Since cmake 3.9, we can use
# find_package(OpenMP). Keep this for older cmake!!!
# ##############################################################################

option(USE_OPENMP "Use OPENMP" OFF)

if(USE_OPENMP)
  if(CMAKE_VERSION VERSION_GREATER 3)
    cmake_policy(SET "CMP0054" NEW)
  endif()

  find_package(OpenMP)
  if(OPENMP_FOUND)
    message(STATUS "INFO: OpenMP found")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  else(OPENMP_FOUND)
    message(WARNING "WARNING: OpenMP was requested but disabled!")
  endif(OPENMP_FOUND)
endif(USE_OPENMP)
