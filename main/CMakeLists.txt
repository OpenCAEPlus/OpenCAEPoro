# Main executable target: testOpenCAEPoro
add_executable(testOpenCAEPoro)
target_sources(testOpenCAEPoro PRIVATE Main.cpp)
target_link_libraries(testOpenCAEPoro PUBLIC
                      OpenCAEPoro
                      ${OPTIONAL_LIBS}
                      fasp
                      ${LAPACK_LIBRARIES}
                      ${BLAS_LIBRARIES} 
                      ${ADD_STDLIBS})
install(TARGETS testOpenCAEPoro
        DESTINATION ${PROJECT_SOURCE_DIR})