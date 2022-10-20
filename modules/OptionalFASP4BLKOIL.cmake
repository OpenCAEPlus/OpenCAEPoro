

##################################################################
# For FASP4BLKOIL
##################################################################

option(USE_FASP4BLKOIL "Use FASP4BLKOIL" OFF)

if(USE_FASP4BLKOIL)

    # set the path to find specific modules
    set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/modules")

    set(FASP4BLKOIL_DIR "${FASP4BLKOIL_DIR}")

    find_package(FASP4BLKOIL)
    if(FASP4BLKOIL_FOUND)
        add_definitions("-DWITH_FASP4BLKOIL=1")
        include_directories(${FASP4BLKOIL_INCLUDE_DIRS})
        set(OPTIONAL_LIBS ${OPTIONAL_LIBS} ${FASP4BLKOIL_LIBRARIES})
    else(FASP4BLKOIL_FOUND)
        message("-- WARNING: FASP4BLKOIL was requested but not supported! Continue without it.")
    endif(FASP4BLKOIL_FOUND)

endif(USE_FASP4BLKOIL)
