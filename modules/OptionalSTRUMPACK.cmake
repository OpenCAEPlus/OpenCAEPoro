##################################################################
# For STRUMPACK
##################################################################

if(USE_STRUMPACK)

    # set the path to find specific modules
    set(STRUMPACK_DIR "${STRUMPACK_DIR}")

    # try to find STRUMPACK
    find_package(STRUMPACK REQUIRED)

    if (STRUMPACK_FOUND)
        add_definitions("-DWITH_STRUMPACK=1")
        include_directories(${STRUMPACK_INCLUDE_DIRS})
        target_link_libraries(${LIBNAME} PUBLIC ${STRUMPACK_LIBRARIES})
    else(STRUMPACK_FOUND)
        message("-- WARNING: STRUMPACK was requested but not supported! Continue without it.")
    endif(STRUMPACK_FOUND)

endif(USE_STRUMPACK)