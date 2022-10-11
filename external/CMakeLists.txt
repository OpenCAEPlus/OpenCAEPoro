include(FetchContent)
# FetchContent_MakeAvailable was not added until CMake 3.14; use our shim
if(${CMAKE_VERSION} VERSION_LESS 3.14)
  include(cmake/add_FetchContent_MakeAvailable.cmake)
endif()

find_package(FASP)

if(FASP_FOUND)
  include_directories(${FASP_INCLUDE_DIRS})
  add_definitions(-D__SOLVER_FASP__)
  add_library(fasp STATIC IMPORTED GLOBAL)
  set_target_properties(
    fasp PROPERTIES IMPORTED_LOCATION ${FASP_LIBRARIES}
                    INTERFACE_INCLUDE_DIRECTORIES ${FASP_INCLUDE_DIRS})
else(FASP_FOUND)
  message("-- INFO FASP was requested but not found!")
  FetchContent_Declare(fasp
                       GIT_REPOSITORY https://github.com/FaspDevTeam/faspsolver.git)
  FetchContent_MakeAvailable(fasp)
endif(FASP_FOUND)

option(BUILD_ECL "Build ECL" OFF)
if(BUILD_ECL)
FetchContent_Declare(
  ecl
  GIT_REPOSITORY https://github.com/billcxx/ecl.git
  GIT_SUBMODULES ""
  GIT_TAG 609fb59)
FetchContent_MakeAvailable(ecl)
endif()
