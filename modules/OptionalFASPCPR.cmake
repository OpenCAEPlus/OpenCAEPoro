# ##############################################################################
# For FASPCPR
# ##############################################################################

option(USE_FASPCPR "Use FASPCPR" OFF)

if(USE_FASPCPR)

  # set the path to find specific modules
  set(FASPCPR_DIR "${FASPCPR_DIR}")

  find_package(FASPCPR)
  if(FASPCPR_FOUND)
    add_library(faspcpr STATIC IMPORTED GLOBAL)
    set_property(
      TARGET faspcpr
      APPEND
      PROPERTY IMPORTED_LOCATION ${FASPCPR_LIBRARIES})
    set_property(
      TARGET faspcpr
      APPEND
      PROPERTY INTERFACE_COMPILE_DEFINITIONS "WITH_FASPCPR=1")
    set_property(
      TARGET faspcpr
      APPEND
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${FASPCPR_INCLUDE_DIRS})
    target_link_libraries(${LIBNAME} PUBLIC faspcpr)
  else(FASPCPR_FOUND)
    message(
      WARNING
        "WARNING: FASPCPR was requested but not found! Continue without it."
    )
  endif(FASPCPR_FOUND)

endif(USE_FASPCPR)
