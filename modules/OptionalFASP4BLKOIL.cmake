# ##############################################################################
# For FASP4BLKOIL
# ##############################################################################

option(USE_FASP4BLKOIL "Use FASP4BLKOIL" OFF)

if(USE_FASP4BLKOIL)

  # set the path to find specific modules
  set(FASP4BLKOIL_DIR "${FASP4BLKOIL_DIR}")

  find_package(FASP4BLKOIL)
  if(FASP4BLKOIL_FOUND)
    add_library(fasp4blkoil STATIC IMPORTED GLOBAL)
    set_property(
      TARGET fasp4blkoil
      APPEND
      PROPERTY IMPORTED_LOCATION ${FASP4BLKOIL_LIBRARIES})
    set_property(
      TARGET fasp4blkoil
      APPEND
      PROPERTY INTERFACE_COMPILE_DEFINITIONS "WITH_FASP4BLKOIL=1")
    set_property(
      TARGET fasp4blkoil
      APPEND
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${FASP4BLKOIL_INCLUDE_DIRS})
    target_link_libraries(${LIBNAME} PUBLIC fasp4blkoil)
  else(FASP4BLKOIL_FOUND)
    message(
      WARNING
        "WARNING: FASP4BLKOIL was requested but not found! Continue without it."
    )
  endif(FASP4BLKOIL_FOUND)

endif(USE_FASP4BLKOIL)
