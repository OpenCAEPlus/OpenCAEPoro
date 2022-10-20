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

  link_directories(${CUDA_DIR}/lib64)
  if(CUDA_FOUND)
    include_directories(${CUDA_INCLUDE_DIRS})
  else(CUDA_FOUND)
    message("-- ERROR: CUDA was requested but not found!")
  endif(CUDA_FOUND)

  # Find FASP4CUDA

  # set the path to find specific modules
  set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/modules")

  set(FASP4CUDA_DIR "${FASP4CUDA_DIR}")

  find_package(FASP4CUDA)
  if(FASP4CUDA_FOUND)
    add_definitions("-DWITH_FASP4CUDA=1")
    include_directories(${FASP4CUDA_INCLUDE_DIRS})
    set(OPTIONAL_LIBS ${OPTIONAL_LIBS} ${FASP4CUDA_LIBRARIES})
  else(FASP4CUDA_FOUND)
    message(
      "-- WARNING: FASP4CUDA was requested but not supported! Continue without it."
    )
  endif(FASP4CUDA_FOUND)

endif(USE_FASP4CUDA)
