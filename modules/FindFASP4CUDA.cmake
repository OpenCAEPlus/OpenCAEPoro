# Once done this will define
#  FASP4CUDA_FOUND        - System has FASP4CUDA
#  FASP4CUDA_DIR          - The FASP4CUDA directory
#  FASP4CUDA_INCLUDE_DIRS - The FASP4CUDA include directories
#  FASP4CUDA_LIBRARIES    - The libraries needed to use FASP4CUDA
#
#  Li Zhao
#  03/20/2022

# message(STATUS "Looking for FASP4CUDA")

set(FASP4CUDA_DIR "${FASP4CUDA_DIR}")

# Check for header file
find_path(FASP4CUDA_INCLUDE_DIRS fasp4cuda.h
   HINTS ${FASP4CUDA_DIR}/include $ENV{FASP4CUDA_DIR}/include ${PROJECT_SOURCE_DIR}/fasp4blkoil/include ${PROJECT_SOURCE_DIR}/fasp/include
   DOC "Directory where the FASP4CUDA header is located")
mark_as_advanced(FASP4CUDA_INCLUDE_DIRS)

# Check for FASP4CUDA library
find_library(FASP4CUDA_LIBRARIES fasp4cuda
    HINTS ${FASP4CUDA_DIR}/lib $ENV{FASP4CUDA_DIR}/lib ${PROJECT_SOURCE_DIR}/fasp4blkoil/lib ${PROJECT_SOURCE_DIR}/fasp/lib
    DOC "The FASP4CUDA library")
mark_as_advanced(FASP4CUDA_LIBRARIES)

set(FASP4CUDA_LIBRARIES ${FASP4CUDA_LIBRARIES})

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FASP4CUDA
    "FASP4CUDA could not be found. Check FASP4CUDA_DIR."
    FASP4CUDA_INCLUDE_DIRS FASP4CUDA_LIBRARIES)

