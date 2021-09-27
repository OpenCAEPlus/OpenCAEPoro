message(STATUS "Checking for dependent packages of 'FASP4BLKOIL'")

set(FASP4BLKOIL_DIR "${FASP4BLKOIL_DIR}")

# Find packages that FASP4BLKOIL depends on
find_package(FASP REQUIRED)

message(STATUS "Looking for FASP4BLKOIL")

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

# Collect libraries
if (FASP_FOUND)
    set(FASP4BLKOIL_LIBRARIES ${FASP4BLKOIL_LIBRARIES})
endif()

# Standard package handling
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FASP4BLKOIL
    "FASP4BLKOIL could not be found. Be sure to set FASP4BLKOIL_DIR."
    FASP4BLKOIL_LIBRARIES FASP4BLKOIL_INCLUDE_DIRS)
