# Main executable target: testOpenCAEPoro
add_executable(testOpenCAEPoro)
target_sources(testOpenCAEPoro PRIVATE Main.cpp)
target_link_libraries(
  testOpenCAEPoro PRIVATE OpenCAEPoro ${OPTIONAL_LIBS} ${FASP_LIBRARIES}
                          ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${ADD_STDLIBS})
install(TARGETS testOpenCAEPoro DESTINATION ${PROJECT_SOURCE_DIR})

if(BUILD_TEST)
  include(CTest)
  add_test(NAME SPE5 COMMAND testOpenCAEPoro ${PROJECT_SOURCE_DIR}/examples/spe5/spe5.data method=IMPEC
                             dtInit=0.1 dtMax=1 dtMin=0.1)
endif()
