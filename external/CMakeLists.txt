include(FetchContent)
# FetchContent_MakeAvailable was not added until CMake 3.14; use our shim
if(${CMAKE_VERSION} VERSION_LESS 3.14)
  include(add_FetchContent_MakeAvailable.cmake)
endif()

# Find Required dependencies 
include(RequiredBLAS)
include(RequiredLAPACK)
include(RequiredFASP)

# OpenCAEPoro library targets

# Find Optional dependencies
include(OptionalOPENMP)
include(OptionalECL)
include(OptionalVTK)
include(OptionalFASP4BLKOIL)
include(OptionalFASP4CUDA)
include(OptionalMUMPS)
include(OptionalPARDISO)
include(OptionalSUPERLU)
include(OptionalUMFPACK)
include(OptionalDoxygen)


# For future, if we have external source to be included, we can put it under
# external and add_subdirectory()
