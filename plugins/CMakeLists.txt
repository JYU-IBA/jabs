#
#    Jaakko's Backscattering Simulator (JaBS)
#    Copyright (C) 2021 - 2023 Jaakko Julin
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#   See LICENSE.txt for the full license.
#

cmake_minimum_required(VERSION 3.15)
project(testplugin
    VERSION 0.0.1
    DESCRIPTION "Minimal test plugin for JABS"
    LANGUAGES C
)
set(CMAKE_C_STANDARD 99)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
configure_file(testplugin_version.h.in testplugin_version.h @ONLY)
find_package(Jibal 0.3.8 REQUIRED)
add_library(testplugin SHARED testplugin.c)
target_link_libraries(testplugin PRIVATE jibal)
install(TARGETS testplugin LIBRARY DESTINATION lib)
