# Once done this will define
#  FASPCPR_FOUND        - System has FASPCPR
#  FASPCPR_DIR          - The FASPCPR directory
#  FASPCPR_INCLUDE_DIRS - The FASPCPR include directories
#  FASPCPR_LIBRARIES    - The libraries needed to use FASPCPR
#
#  Chensong Zhang, Li Zhao
#  12/11/2022

# message(STATUS "Looking for FASPCPR")

set(FASPCPR_DIR "${FASPCPR_DIR}")

# Check for header file
find_path(FASPCPR_INCLUDE_DIRS faspcpr.h
   HINTS ${FASPCPR_DIR}/include $ENV{FASPCPR_DIR}/include ${PROJECT_SOURCE_DIR}/fasp/include
   DOC "Directory where the FASPCPR header is located")
mark_as_advanced(FASPCPR_INCLUDE_DIRS)

# Check for FASPCPR library
find_library(FASPCPR_LIBRARIES faspcpr
    HINTS ${FASPCPR_DIR}/lib $ENV{FASPCPR_DIR}/lib ${PROJECT_SOURCE_DIR}/fasp/lib
    DOC "The FASPCPR library")
mark_as_advanced(FASPCPR_LIBRARIES)

set(FASPCPR_LIBRARIES ${FASPCPR_LIBRARIES})

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FASPCPR
    "FASPCPR could not be found. Check FASPCPR_DIR."
    FASPCPR_INCLUDE_DIRS FASPCPR_LIBRARIES)
