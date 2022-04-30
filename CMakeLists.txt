#
#    Jaakko's Backscattering Simulator (JaBS)
#    Copyright (C) 2021 Jaakko Julin
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#   See LICENSE.txt for the full license.
#
cmake_minimum_required(VERSION 3.15)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
if(WIN32 AND DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
      CACHE STRING "")
  message(STATUS "CMake toolchain file ${CMAKE_TOOLCHAIN_FILE}")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release)
endif()

add_subdirectory(gitwatcher)

#Read version from file
file(READ "${CMAKE_CURRENT_LIST_DIR}/version.txt" BUILD_VERSION)
string(STRIP "${BUILD_VERSION}" BUILD_VERSION)


project(jabs
        VERSION "${BUILD_VERSION}"
        DESCRIPTION "Jaakko's Backscattering Simulator (JaBS)"
        LANGUAGES C
        )
set(CMAKE_C_STANDARD 99)

add_subdirectory(src)
