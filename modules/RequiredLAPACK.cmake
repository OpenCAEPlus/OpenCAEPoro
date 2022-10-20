find_package(LAPACK REQUIRED)

# include_directories(${LAPACK_INCLUDE_DIRS}) message(STATUS 1231231
# ${LAPACK_LIBRARIES})
add_library(lapack INTERFACE IMPORTED GLOBAL)
set_property(
  TARGET lapack
  APPEND
  PROPERTY INTERFACE_LINK_LIBRARIES ${LAPACK_LIBRARIES})
# INTERFACE_INCLUDE_DIRECTORIES ${LAPACK_INCLUDE_DIRS}

target_link_libraries(${LIBNAME} PUBLIC lapack)
