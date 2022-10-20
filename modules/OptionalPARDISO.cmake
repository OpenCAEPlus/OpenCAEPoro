##################################################################
# For Intel MKL PARDISO
##################################################################

option(USE_PARDISO "Use PARDISO" OFF)

if(USE_PARDISO)

    # set the path to find specific modules
    set(MKL_DIR "${MKL_DIR}")

    # try to find MKL
    find_package(MKL)

    if (MKL_FOUND)
        add_definitions("-DWITH_PARDISO=1")
        include_directories(${MKL_INCLUDE_DIRS})
        target_link_libraries(${LIBNAME} PUBLIC ${MKL_LIBRARIES})
    else(MKL_FOUND)
        message("-- WARNING: Intel MKL was requested but not supported! Continue without it.")
    endif(MKL_FOUND)

endif(USE_PARDISO)
