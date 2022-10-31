# ##############################################################################
# For ECL
# ##############################################################################

option(USE_ECL "Use ECL" OFF)

if(USE_ECL)
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL Intel OR ${CMAKE_CXX_COMPILER_ID}
                                                STREQUAL IntelLLVM)
    find_package(ECL)
    message(STATUS "For intel compiler, you have to install ECL by yourself")
  else()
    FetchContent_Declare(
      ECL
      GIT_REPOSITORY https://github.com/billcxx/ecl.git
      GIT_SUBMODULES ""
      GIT_TAG 609fb59)
    FetchContent_MakeAvailable(ECL)
  endif()

  target_link_libraries(${LIBNAME} PUBLIC ECL)

endif()
