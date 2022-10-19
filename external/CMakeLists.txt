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
  FetchContent_Declare(
    fasp GIT_REPOSITORY https://github.com/FaspDevTeam/faspsolver.git)
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

set(LINK_VTK_COMPONENTS CommonCore IOXML)
find_package(VTK COMPONENTS ${LINK_VTK_COMPONENTS})
if(VTK_FOUND)
  message(STATUS "Found vtk in system")
else()
  message(STATUS "Going to use vtk in local folder")
  FetchContent_Declare(
    vtk
    GIT_REPOSITORY https://gitlab.kitware.com/vtk/vtk.git
    GIT_TAG v9.2.2)
  FetchContent_GetProperties(vtk)
  if(NOT vtk_POPULATED)
    message(STATUS "Auto download VTK")
    # Fetch the content using previously declared details
    FetchContent_Populate(vtk)
  endif()
  message(STATUS "Build and install VTK")
  add_custom_target(
    vtk ALL
    COMMAND
      ${CMAKE_COMMAND} -S "${vtk_SOURCE_DIR}" -B . -G "${CMAKE_GENERATOR}"
      -DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
      -DCMAKE_CXX_COMPILER="${CMAKE_CXX_COMPILER}"
      -DCMAKE_C_COMPILER="${CMAKE_C_COMPILER}"
      -DCMAKE_C_FLAGS="-w"
      -DCMAKE_CXX_FLAGS="-w"
      -DCMAKE_INSTALL_PREFIX="${PROJECT_SOURCE_DIR}/out/vtk/${CMAKE_SYSTEM_NAME}-${CMAKE_CXX_COMPILER_ID}-${CMAKE_BUILD_TYPE}"
    COMMAND ${CMAKE_COMMAND} --build . -j
    COMMAND ${CMAKE_COMMAND} --install .
    WORKING_DIRECTORY ${vtk_BINARY_DIR}
    USES_TERMINAL)

  message(
    STATUS "VTK install path "
           ${CMAKE_SYSTEM_NAME}-${CMAKE_CXX_COMPILER_ID}-${CMAKE_BUILD_TYPE})
  set(VTK_DIR
      "${PROJECT_SOURCE_DIR}/out/vtk/${CMAKE_SYSTEM_NAME}-${CMAKE_CXX_COMPILER_ID}-${CMAKE_BUILD_TYPE}/lib/cmake/vtk-9.2"
  )
  find_package(VTK COMPONENTS ${LINK_VTK_COMPONENTS})
  message(STATUS ${VTK_FOUND} ${VTK_LIBRARIES})
endif()
