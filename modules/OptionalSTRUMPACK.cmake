# ##############################################################################
# For STRUMPACK
# ##############################################################################

if(USE_STRUMPACK)

  # set the path to find specific modules
  set(STRUMPACK_DIR "${STRUMPACK_DIR}")

  # try to find STRUMPACK
  find_package(STRUMPACK REQUIRED)

  if(STRUMPACK_FOUND)
    message(STATUS "INFO: STRUMPACK found")
    add_library(strumpack INTERFACE IMPORTED GLOBAL)
    set_property(
      TARGET strumpack
      APPEND
      PROPERTY INTERFACE_LINK_LIBRARIES ${STRUMPACK_LIBRARIES})
    set_property(
      TARGET strumpack
      APPEND
      PROPERTY INTERFACE_COMPILE_DEFINITIONS "WITH_STRUMPACK=1")
    set_property(
      TARGET strumpack
      APPEND
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${STRUMPACK_INCLUDE_DIRS})
    target_link_libraries(${LIBNAME} PUBLIC strumpack)
  else(STRUMPACK_FOUND)
    message(
      WARNING
        "WARNING: STRUMPACK was requested but not found! Continue without it."
    )
  endif(STRUMPACK_FOUND)

endif(USE_STRUMPACK)
