# ##############################################################################
# For FASP4BLKOIL
# ##############################################################################

option(USE_FASP4BLKOIL "Use FASP4BLKOIL" OFF)

if(USE_FASP4BLKOIL)

  # set the path to find specific modules
  set(FASP4BLKOIL_DIR "${FASP4BLKOIL_DIR}")

  find_package(FASP4BLKOIL)
  if(FASP4BLKOIL_FOUND)
    message(STATUS "INFO: FASP4BLKOIL found")
    add_library(fasp4blkoil STATIC IMPORTED GLOBAL)
    set_property(
      TARGET fasp4blkoil
      APPEND
      PROPERTY IMPORTED_LOCATION ${FASP4BLKOIL_LIBRARIES})
    set_property(TARGET fasp4blkoil APPEND COMPILE_DEFINITIONS
                                "-DWITH_FASP4BLKOIL=1")
    set_property(TARGET fasp4blkoil APPEND INTERFACE_INCLUDE_DIRECTORIES
                                ${FASP4BLKOIL_INCLUDE_DIRS})
    target_link_libraries(${LIBNAME} PUBLIC fasp4blkoil)

    # add_definitions("-DWITH_FASP4BLKOIL=1")
    # include_directories(${FASP4BLKOIL_INCLUDE_DIRS}) set(OPTIONAL_LIBS
    # ${OPTIONAL_LIBS} ${FASP4BLKOIL_LIBRARIES})
  else(FASP4BLKOIL_FOUND)
    message(
      WARNING
        "WARNING: FASP4BLKOIL was requested but not supported! Continue without it."
    )
  endif(FASP4BLKOIL_FOUND)

  target_link_libraries(${LIBNAME} PUBLIC ${FASP4BLKOIL_LIBRARIES})

endif(USE_FASP4BLKOIL)
