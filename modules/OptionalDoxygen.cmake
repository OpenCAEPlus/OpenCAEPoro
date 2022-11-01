

##################################################################
# For Doxygen
##################################################################
option(USE_DOXYGEN "Use DOXYGEN" OFF)

if(USE_DOXYGEN)

    find_package(Doxygen)

    if(DOXYGEN_FOUND)
        if(EXISTS ${FASP_SOURCE_DIR}/doc/fasp.Doxygen.cnf.in)
            configure_file(
                    ${FASP_SOURCE_DIR}/doc/fasp.Doxygen.cnf.in
                    ${CMAKE_CURRENT_BINARY_DIR}/fasp.Doxygen.cnf @ONLY)
            set(DOXY_EXEC "${DOXYGEN_EXECUTABLE}")
            if(DOXYWIZARD)
                find_program(WIZARD doxywizard)
                if(APPLE AND (NOT WIZARD))
                    find_program(WIZARD
                            /Applications/Doxygen.app/Contents/MacOS/Doxywizard)
                endif()
                if(WIZARD)
                    set(DOXY_EXEC "${WIZARD}")
                endif()
            endif(DOXYWIZARD)
            add_custom_target(docs ${DOXY_EXEC}
                    ${CMAKE_CURRENT_BINARY_DIR}/fasp.Doxygen.cnf
                    WORKING_DIRECTORY
                    "${CMAKE_CURRENT_BINARY_DIR}"
                    COMMENT
                    "Generating FASP documentation (Doxygen)"
                    VERBATIM)
        else(EXISTS ${FASP_SOURCE_DIR}/doc/fasp.Doxygen.cnf.in)
            message(
                WARNING
                "-- WARNING: Doxygen configuration file cannot be found! Continue without it.")
        endif(EXISTS ${FASP_SOURCE_DIR}/doc/fasp.Doxygen.cnf.in)
    endif(DOXYGEN_FOUND)

endif(USE_DOXYGEN)
