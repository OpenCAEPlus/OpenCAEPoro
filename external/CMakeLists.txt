include(FetchContent)
# FetchContent_MakeAvailable was not added until CMake 3.14; use our shim
if(${CMAKE_VERSION} VERSION_LESS 3.14)
  include(add_FetchContent_MakeAvailable.cmake)
endif()

# Find required dependencies 
include(RequiredBLAS)
include(RequiredLAPACK)
include(RequiredFASP)

# Find optional dependencies
include(OptionalOPENMP)
include(OptionalFASP4BLKOIL)
include(OptionalFASP4CUDA)
include(OptionalMUMPS)
include(OptionalPARDISO)
include(OptionalSUPERLU)
include(OptionalUMFPACK)
include(OptionalDoxygen)
include(OptionalECL)
include(OptionalVTK)

# For future, if we have external source to be included, we can put it under
# external and add_subdirectory()
