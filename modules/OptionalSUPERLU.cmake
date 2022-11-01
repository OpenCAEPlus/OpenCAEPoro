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
    set_property(
      TARGET superlu
      APPEND
      PROPERTY INTERFACE_COMPILE_DEFINITIONS "WITH_SuperLU=1")
    set_property(
      TARGET superlu
      APPEND
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${SUPERLU_INCLUDE_DIRS})
    target_link_libraries(${LIBNAME} PUBLIC superlu)
  else(SUPERLU_FOUND)
    message(
      WARNING
        "WARNING: SuperLU was requested but not found! Continue without it."
    )
  endif(SUPERLU_FOUND)

endif(USE_SUPERLU)
