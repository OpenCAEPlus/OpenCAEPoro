option(USE_ECL "Use ECL" OFF)

if(USE_ECL)
  FetchContent_Declare(
    ecl
    GIT_REPOSITORY https://github.com/billcxx/ecl.git
    GIT_SUBMODULES ""
    GIT_TAG 609fb59)
  FetchContent_MakeAvailable(ecl)

  target_link_libraries(${LIBNAME} PUBLIC ecl)

endif()
