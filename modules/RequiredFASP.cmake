# ##############################################################################
# For FASPSOLVER
# ##############################################################################

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

target_link_libraries(${LIBNAME} PUBLIC fasp)
