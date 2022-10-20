##################################################################
# For UMFPACK
##################################################################

option(USE_UMFPACK "Use UMFPACK" OFF)

if(USE_UMFPACK)

    # set some path to the UMFPACK pacakge
    # metis is not part of suitesparse, so theremay be also some other metis dir.
    set(METIS_DIR "${SUITESPARSE_DIR}")

    find_package(UMFPACK)
    if (UMFPACK_FOUND)
        add_definitions("-DWITH_UMFPACK=1")
        include_directories(${UMFPACK_INCLUDE_DIRS})
        target_link_libraries(${LIBNAME} PUBLIC ${UMFPACK_LIBRARIES})
    else(UMFPACK_FOUND)
        message("-- WARNING: UMFPACK was requested but not supported! Continue without it.")
    endif(UMFPACK_FOUND)

endif(USE_UMFPACK)
