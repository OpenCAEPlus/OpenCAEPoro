# ##############################################################################
# For UMFPACK
# ##############################################################################

option(USE_UMFPACK "Use UMFPACK" OFF)

if(USE_UMFPACK)

  # set some path to the UMFPACK pacakge metis is not part of suitesparse, so
  # theremay be also some other metis dir.
  set(SUITESPARSE_DIR "${SUITESPARSE_DIR}")
  find_package(UMFPACK)
  if(UMFPACK_FOUND)
    message(STATUS "INFO: UMFPACK found")
    add_library(umfpack INTERFACE IMPORTED GLOBAL)
    set_property(
      TARGET umfpack
      APPEND
      PROPERTY INTERFACE_LINK_LIBRARIES ${UMFPACK_LIBRARIES})
    set_property(
      TARGET umfpack
      APPEND
      PROPERTY INTERFACE_COMPILE_DEFINITIONS "WITH_UMFPACK=1")
    set_property(
      TARGET umfpack
      APPEND
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${UMFPACK_INCLUDE_DIRS})
    target_link_libraries(${LIBNAME} PUBLIC umfpack)
  else(UMFPACK_FOUND)
    message(
      WARNING
        "WARNING: UMFPACK was requested but not found! Continue without it."
    )
  endif(UMFPACK_FOUND)

endif(USE_UMFPACK)
