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
include(version.cmake)
project(jabs
        VERSION "${BUILD_VERSION}"
        DESCRIPTION "Jaakko's Backscattering Simulator (JaBS)"
        LANGUAGES C
        )

if(WIN32 AND DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
      CACHE STRING "")
  message(STATUS "CMake toolchain file ${CMAKE_TOOLCHAIN_FILE}")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release)
endif()


option(DEBUG_MODE "Debug mode" OFF)
if(DEBUG_MODE)
    if (MSVC)
        add_compile_options(/W4 /WX /DDEBUG)
    else()
        add_compile_options(-Wall -Wextra -pedantic -DDEBUG)
    endif()
endif()

option(PLUGINS "Support for plugins" OFF)
if(PLUGINS)
    if (MSVC)
        #No Windows support in the first phase
    else()
        add_compile_options(-DJABS_PLUGIN_SUPPORT)
    endif()
endif()

add_subdirectory(gitwatcher)

add_subdirectory(src)

add_subdirectory(plugins)
