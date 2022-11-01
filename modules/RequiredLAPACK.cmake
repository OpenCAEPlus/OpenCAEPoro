# ##############################################################################
# For LAPACK
# ##############################################################################

find_package(LAPACK REQUIRED)

add_library(lapack INTERFACE IMPORTED GLOBAL)
set_property(
  TARGET lapack
  APPEND
  PROPERTY INTERFACE_LINK_LIBRARIES ${LAPACK_LIBRARIES})

target_link_libraries(${LIBNAME} PUBLIC lapack)
