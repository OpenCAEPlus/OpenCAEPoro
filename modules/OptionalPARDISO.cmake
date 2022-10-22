# ##############################################################################
# For Intel MKL PARDISO
# ##############################################################################

option(USE_PARDISO "Use PARDISO" OFF)

if(USE_PARDISO)

  # set the path to find specific modules
  set(MKL_DIR "${MKL_DIR}")

  # try to find MKL
  find_package(MKL)

  if(MKL_FOUND)
    message(STATUS "INFO: Intel MKL found")
    add_library(pardiso INTERFACE IMPORTED GLOBAL)
    set_property(TARGET pardiso APPEND PROPERTY INTERFACE_LINK_LIBRARIES
                                ${MKL_LIBRARIES})
    set_property(TARGET pardiso APPEND PROPERTY COMPILE_DEFINITIONS "-DWITH_PARDISO=1")
    set_property(TARGET pardiso APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                ${MKL_INCLUDE_DIRS})

    target_link_libraries(${LIBNAME} PUBLIC pardiso)

    # add_definitions("-DWITH_PARDISO=1")
    # include_directories(${MKL_INCLUDE_DIRS})
  else(MKL_FOUND)
    message(
      WARNING "WARNING: Intel MKL was requested but not supported! Continue without it."
    )
  endif(MKL_FOUND)

endif(USE_PARDISO)
