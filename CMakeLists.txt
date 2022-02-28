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
if(WIN32 AND DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
      CACHE STRING "")
  message(STATUS "CMake toolchain file ${CMAKE_TOOLCHAIN_FILE}")
endif()
project(jabs
        VERSION 0.5.2
        DESCRIPTION "Jaakko's Backscattering Simulator (JaBS)"
        LANGUAGES C
        )

set(CMAKE_C_STANDARD 99)

configure_file(version.h.in version.h @ONLY)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(GSL 2.6 REQUIRED)
find_package(Jibal 0.3.5 REQUIRED)

if(WIN32)
    find_path(GETOPT_INCLUDE_DIR getopt.h)
    find_library(GETOPT_LIBRARY getopt)
endif()

option(DEBUG_MODE_ENABLE "Enable debug mode" OFF)
if(DEBUG_MODE_ENABLE)
    if (MSVC)
        add_compile_options(/W4 /WX /DDEBUG)
    else()
        add_compile_options(-Wall -Wextra -pedantic -DDEBUG)
    endif()
endif()

add_executable(jabs main.c sample.c sample.h brick.c ion.c ion.h simulation.c simulation.h reaction.c reaction.h
        options.c options.h spectrum.c spectrum.h fit.c fit.h rotate.c rotate.h detector.c detector.h jabs.c jabs.h
        roughness.c roughness.h script.c script.h generic.c generic.h message.c message.h aperture.c aperture.h
        geostragg.c geostragg.h script_command.c script_command.h script_session.c script_session.h script_file.c
        script_file.h calibration.c calibration.h script_generic.h)

target_include_directories(jabs PRIVATE ${GETOPT_INCLUDE_DIR})

target_link_libraries(jabs
    PRIVATE jibal
    PRIVATE GSL::gsl
    PRIVATE "$<$<BOOL:${UNIX}>:m>"
    ${GETOPT_LIBRARY}
    )

install(TARGETS jabs RUNTIME DESTINATION bin)
