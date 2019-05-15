cmake_minimum_required(VERSION 3.1)

project(
  "Digital Signal Analysis C library"
  VERSION 1.0
  LANGUAGES C
)

if(NOT DEFINED CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE "RelWithDebInfo")
endif()

option(DSAC_DEBUG  "Activate debug build" OFF)
option(DSAC_OPENMP "Activate the OpenMP parallelization" OFF)

# CFLAGS
# assume it's a GNU-like compiler (GCC, clang)
set(CMAKE_C_FLAGS                 "-Wall -Wextra -pedantic")
set(CMAKE_C_FLAGS_DEBUG           "-O0 -g")
set(CMAKE_C_FLAGS_MINSIZEREL      "-Os -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE         "-O3 -DNDEBUG")
set(CMAKE_C_FLAGS_RELWITHDEBINFO  "-O2 -g")

if (DSAC_DEBUG)
  message(STATUS "Debug enable")
  set(CMAKE_BUILD_TYPE Debug)
endif()

# Option for OpenMP
if (DSAC_OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -DDSA_OPENMP")
    message(STATUS "OpenMP enable")
  else()
    message(SEND_ERROR "OpenMP not found!")
  endif()
else()
  message(STATUS "OpenMP disable")
endif()

# Depends
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules")
find_package(FFTW3 REQUIRED)

# Include library sources
add_subdirectory("src")

# Include GNUInstallDirs (if Unix system)
if(UNIX)
  include(GNUInstallDirs)
endif()

# Rules for install
install(TARGETS dsac
  EXPORT dsacConfig
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/dsac
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING
  PATTERN "*.h"
  PATTERN ".swp" EXCLUDE
)

install(
  EXPORT dsacConfig
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/dsac"
)
