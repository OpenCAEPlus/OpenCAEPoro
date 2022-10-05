# Once done this will define
#  FASP_FOUND        - System has FASP
#  FASP_DIR          - The FASP directory
#  FASP_INCLUDE_DIRS - The FASP include directories
#  FASP_LIBRARIES    - The libraries needed to use FASP
#
#  Chensong Zhang
#  01/18/2022

set(FASP_DIR "${FASP_DIR}")

# Check for header file
find_path(FASP_INCLUDE_DIRS fasp.h
   HINTS ${FASP_DIR}/include $ENV{FASP_DIR}/include ${PROJECT_SOURCE_DIR}/fasp/include
   DOC "Directory where the FASP header is located")
mark_as_advanced(FASP_INCLUDE_DIRS)

# Check for FASP library
find_library(FASP_LIBRARIES fasp
    HINTS ${FASP_DIR}/lib $ENV{FASP_DIR}/lib ${PROJECT_SOURCE_DIR}/fasp/lib
    DOC "The FASP library")
mark_as_advanced(FASP_LIBRARIES)

# Collect libraries
set(FASP_LIBRARIES ${FASP_LIBRARIES})

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FASP
    "FASP could not be found. Check FASP_DIR."
    FASP_LIBRARIES FASP_INCLUDE_DIRS)
