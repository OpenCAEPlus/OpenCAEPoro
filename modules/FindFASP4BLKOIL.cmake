# Once done this will define
#  FASP4BLKOIL_FOUND        - System has FASP4BLKOIL
#  FASP4BLKOIL_DIR          - The FASP4BLKOIL directory
#  FASP4BLKOIL_INCLUDE_DIRS - The FASP4BLKOIL include directories
#  FASP4BLKOIL_LIBRARIES    - The libraries needed to use FASP4BLKOIL
#
#  Chensong Zhang
#  01/18/2022

# message(STATUS "Looking for FASP4BLKOIL")

set(FASP4BLKOIL_DIR "${FASP4BLKOIL_DIR}")

# Check for header file
find_path(FASP4BLKOIL_INCLUDE_DIRS fasp4blkoil.h
   HINTS ${FASP4BLKOIL_DIR}/include $ENV{FASP4BLKOIL_DIR}/include ${PROJECT_SOURCE_DIR}/fasp/include
   DOC "Directory where the FASP4BLKOIL header is located")
mark_as_advanced(FASP4BLKOIL_INCLUDE_DIRS)

# Check for FASP4BLKOIL library
find_library(FASP4BLKOIL_LIBRARIES fasp4blkoil
    HINTS ${FASP4BLKOIL_DIR}/lib $ENV{FASP4BLKOIL_DIR}/lib ${PROJECT_SOURCE_DIR}/fasp/lib
    DOC "The FASP4BLKOIL library")
mark_as_advanced(FASP4BLKOIL_LIBRARIES)

set(FASP4BLKOIL_LIBRARIES ${FASP4BLKOIL_LIBRARIES})

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FASP4BLKOIL
    "FASP4BLKOIL could not be found. Check FASP4BLKOIL_DIR."
    FASP4BLKOIL_INCLUDE_DIRS FASP4BLKOIL_LIBRARIES)
