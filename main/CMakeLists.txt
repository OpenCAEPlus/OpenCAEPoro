# Main executable target: testOpenCAEPoro
add_executable(testOpenCAEPoro)
target_sources(testOpenCAEPoro PRIVATE Main.cpp)
target_link_libraries(testOpenCAEPoro PUBLIC OpenCAEPoro ${ADD_STDLIBS})
install(TARGETS testOpenCAEPoro DESTINATION ${PROJECT_SOURCE_DIR})

if(BUILD_TEST)
  include(CTest)
  add_test(
    NAME SPE5
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples/spe5/
    COMMAND testOpenCAEPoro spe5.data method=IMPEC dtInit=0.1 dtMax=1 dtMin=0.1)

  add_test(
    NAME SPE1A
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/examples/spe1a/
    COMMAND testOpenCAEPoro spe1a.data method=IMPEC dtInit=0.1 dtMax=1
            dtMin=0.1)
endif()
