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
project(jabs
        VERSION 0.1.4
        DESCRIPTION "Jaakko's Backscattering Simulator (JaBS)"
        LANGUAGES C
        )

set(CMAKE_C_STANDARD 99)

configure_file(version.h.in version.h @ONLY)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(GSL 2.6 REQUIRED)
find_package(Jibal 0.2.6 REQUIRED)

if(WIN32)
    find_path(GETOPT_INCLUDE_DIR getopt.h)
    find_library(GETOPT_LIBRARY getopt)
endif()

option(DEBUG_OUTPUT_ENABLE "Enable debug output" OFF)
if(DEBUG_OUTPUT_ENABLE)
    add_compile_options(-DDEBUG)
endif()

add_executable(jabs main.c sample.c sample.h brick.c ion.c ion.h simulation.c simulation.h
        reaction.c reaction.h layers.c layers.h options.c options.h spectrum.c spectrum.h fit.c fit.h)
target_include_directories(jabs PRIVATE ${GETOPT_INCLUDE_DIR})

target_link_libraries(jabs
    PRIVATE Jibal 
    PRIVATE GSL::gsl
    PRIVATE "$<$<BOOL:${UNIX}>:m>" ${GETOPT_LIBRARY}
    )

install(TARGETS jabs RUNTIME DESTINATION bin)
