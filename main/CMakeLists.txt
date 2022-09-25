# Main executable target: testOpenCAEPoro
add_executable(testOpenCAEPoro)
target_sources(testOpenCAEPoro PRIVATE Main.cpp)
target_link_libraries(testOpenCAEPoro PRIVATE
                      OpenCAEPoro
                      ${OPTIONAL_LIBS}
                      ${FASP_LIBRARIES}
                      ${LAPACK_LIBRARIES}
                      ${BLAS_LIBRARIES} 
                      ${ADD_STDLIBS})
install(TARGETS testOpenCAEPoro
        DESTINATION ${PROJECT_SOURCE_DIR})