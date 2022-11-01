# ##############################################################################
# For BLAS
# ##############################################################################

find_package(BLAS REQUIRED)

# include_directories(${BLAS_INCLUDE_DIRS})
add_library(blas INTERFACE IMPORTED GLOBAL)
set_property(
  TARGET blas
  APPEND
  PROPERTY INTERFACE_LINK_LIBRARIES ${BLAS_LIBRARIES})
# INTERFACE_INCLUDE_DIRECTORIES ${BLAS_INCLUDE_DIRS}

target_link_libraries(${LIBNAME} PUBLIC blas)
