
##################################################################
# For MUMPS
##################################################################
option(USE_MUMPS "Use MUMPS" OFF)

if(USE_MUMPS)

    # set the path to find specific modules
    set(MUMPS_DIR "${MUMPS_DIR}")

    # try to find MUMPS and METIS (as dependency)
    find_package(METIS)
    find_package(MUMPS)

    if (MUMPS_FOUND)
        add_definitions("-DWITH_MUMPS=1")
        include_directories(${MUMPS_INCLUDE_DIRS})
    else(MUMPS_FOUND)
        message("-- WARNING: MUMPS was requested but not supported! Continue without it.")
    endif(MUMPS_FOUND)

endif(USE_MUMPS)