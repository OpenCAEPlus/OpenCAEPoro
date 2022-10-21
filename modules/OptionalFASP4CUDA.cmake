# Set optional libraries in OPTIONAL_LIBS

# ##############################################################################
# For FASP4CUDA
# ##############################################################################

option(USE_FASP4CUDA "Use FASP4CUDA" OFF)

# Find CUDA
if(USE_FASP4CUDA)

  set(CUDA_DIR /usr/local/cuda-10.2)
  set(CUDA_ROOT /usr/local/cuda-10.2) # Set the root directory of CUDA 10.2
                                      # version
  find_package(CUDA REQUIRED)
  # link_directories(${CUDA_DIR}/lib64)
  if(CUDA_FOUND)
    # include_directories(${CUDA_INCLUDE_DIRS})
    message(STATUS "INFO: CUDA found")
    add_library(cuda INTERFACE IMPORTED GLOBAL)
    set_property(TARGET cuda APPEND LINK_DIRECTORIES ${CUDA_DIR}/lib64)
    set_property(TARGET cuda APPEND INTERFACE_LINK_LIBRARIES cublas cusparse
                             cudart cudadevrt)
    set_property(TARGET cuda APPEND INTERFACE_INCLUDE_DIRECTORIES
                             ${CUDA_INCLUDE_DIRS})
  else(CUDA_FOUND)
    message(FATAL_ERROR "ERROR: CUDA was requested but not found!")
  endif(CUDA_FOUND)

  # Find FASP4CUDA

  # set the path to find specific modules
  set(FASP4CUDA_DIR "${FASP4CUDA_DIR}")

  find_package(FASP4CUDA)
  if(FASP4CUDA_FOUND)
    message(STATUS "INFO: FASP4CUDA found")
    add_library(fasp4cuda INTERFACE IMPORTED GLOBAL)
    set_property(TARGET fasp4cuda APPEND INTERFACE_LINK_LIBRARIES cuda
                                  ${FASP4CUDA_LIBRARIES})
    set_property(TARGET fasp4cuda APPEND COMPILE_DEFINITIONS
                                  "-DWITH_FASP4CUDA=1")
    set_property(TARGET fasp4cuda APPEND INTERFACE_INCLUDE_DIRECTORIES
                                  ${FASP4CUDA_INCLUDE_DIRS})
    # add_definitions("-DWITH_FASP4CUDA=1")
    # include_directories(${FASP4CUDA_INCLUDE_DIRS}) set(OPTIONAL_LIBS
    # ${OPTIONAL_LIBS} ${FASP4CUDA_LIBRARIES})
    target_link_libraries(${LIBNAME} PUBLIC fasp4cuda)
  else(FASP4CUDA_FOUND)
    message(
      WANRING
      "WARNING: FASP4CUDA was requested but not supported! Continue without it."
    )
  endif(FASP4CUDA_FOUND)

endif(USE_FASP4CUDA)
