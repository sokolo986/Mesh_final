# Install script for directory: /home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/src/mlpack/methods

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "1")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/cf/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/det/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/emst/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/fastmks/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/gmm/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/hmm/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/kernel_pca/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/kmeans/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/lars/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/linear_regression/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/local_coordinate_coding/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/logistic_regression/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/lsh/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/naive_bayes/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/nca/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/neighbor_search/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/nmf/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/pca/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/radical/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/range_search/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/rann/cmake_install.cmake")
  INCLUDE("/home/sokolo/Desktop/CS207/final/Mesh_final/mlpack-1.0.8/build/src/mlpack/methods/sparse_coding/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

