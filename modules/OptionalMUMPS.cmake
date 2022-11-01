# ##############################################################################
# For MUMPS
# ##############################################################################

option(USE_MUMPS "Use MUMPS" OFF)

if(USE_MUMPS)

  # set the path to find specific modules
  set(MUMPS_DIR "${MUMPS_DIR}")

  # try to find MUMPS and METIS (as dependency)
  find_package(METIS)
  find_package(MUMPS)

  if(MUMPS_FOUND)
    message(STATUs "INFO: MUMPS found")
    add_library(mumps INTERFACE IMPORTED GLOBAL)
    set_property(
      TARGET mumps
      APPEND
      PROPERTY INTERFACE_LINK_LIBRARIES ${METIS_LIBRARIES} ${MUMPS_LIBRARIES})
    set_property(
      TARGET mumps
      APPEND
      PROPERTY INTERFACE_COMPILE_DEFINITIONS "WITH_MUMPS=1")
    set_property(
      TARGET mumps
      APPEND
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${MUMPS_INCLUDE_DIRS})
    target_link_libraries(${LIBNAME} PUBLIC mumps)
  else(MUMPS_FOUND)
    message(
      WARNING
        "WARNING: MUMPS was requested but not found! Continue without it.")
  endif(MUMPS_FOUND)

endif(USE_MUMPS)
