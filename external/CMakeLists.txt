include(FetchContent)
# FetchContent_MakeAvailable was not added until CMake 3.14; use our shim
if(${CMAKE_VERSION} VERSION_LESS 3.14)
  include(add_FetchContent_MakeAvailable.cmake)
endif()

# Find BLAS and LAPACK implementations
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Find FASP, will be installed if not found
include(RequiredFASP)

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
