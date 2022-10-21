# ##############################################################################
# For SuperLU
# ##############################################################################
option(USE_SUPERLU "Use SUPERLU" OFF)

if(USE_SUPERLU)

  # set the path to find specific modules
  set(SUPERLU_DIR "${SUPERLU_DIR}")

  # try to find SuperLU
  find_package(SUPERLU)

  if(SUPERLU_FOUND)
  message(STATUS "INFO: SuperLU found")
    add_library(superlu INTERFACE IMPORTED GLOBAL)
    set_property(TARGET superlu APPEND PROPERTY COMPILE_DEFINITIONS "-DWITH_SuperLU=1")
    set_property(TARGET superlu APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                ${SUPERLU_INCLUDE_DIRS})
    target_link_libraries(${LIBNAME} PUBLIC superlu)

    # add_definitions("-DWITH_SuperLU=1")
    # include_directories(${SUPERLU_INCLUDE_DIRS})
  else(SUPERLU_FOUND)
    message(
      WARNING
        "WARNING: SuperLU was requested but not supported! Continue without it."
    )
  endif(SUPERLU_FOUND)

endif(USE_SUPERLU)
