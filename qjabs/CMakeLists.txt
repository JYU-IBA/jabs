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

cmake_minimum_required(VERSION 3.19)
if(WIN32 AND DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
            CACHE STRING "")
    message(STATUS "CMake toolchain file ${CMAKE_TOOLCHAIN_FILE}")
endif()
include(../version.cmake)
project(qjabs
    VERSION "${BUILD_VERSION}"
    DESCRIPTION "Jaakko's Backscattering Simulator (JaBS) for Qt"
    LANGUAGES C CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(../src)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MACOSX_BUNDLE TRUE)

include(../gitwatcher/git_watcher_config.cmake)

configure_file(version.h.in version.h @ONLY)

find_package(Qt6 COMPONENTS Widgets REQUIRED PrintSupport REQUIRED Svg REQUIRED)

if(${QT_VERSION_MAJOR} NOT EQUAL 6)
    message("Qt major version is ${QT_VERSION_MAJOR}. This program is only for Qt6.")
endif()

find_package(GSL 2.6 REQUIRED)
find_package(Jibal 0.3.8 REQUIRED)
find_package(LibXml2 REQUIRED)
find_package(OpenMP)


if(WIN32)
    find_path(GETOPT_INCLUDE_DIR getopt.h)
    find_library(GETOPT_LIBRARY getopt)
endif()

option(JABS_PLUGINS "Support for plugins" OFF)
if(JABS_PLUGINS)
    if (MSVC)
        #No Windows support in the first phase
    else()
        add_compile_options(-DJABS_PLUGINS)
    endif()
endif()

set(PROJECT_SOURCES
        ../src/sample.c
        ../src/brick.c
        ../src/ion.c
        ../src/simulation.c
        ../src/reaction.c
        ../src/sample.c
        ../src/spectrum.c
        ../src/fit.c
        ../src/rotate.c
        ../src/jabs.c
        ../src/generic.c
        ../src/detector.c
        ../src/roughness.c
        ../src/script.c
        ../src/jabs.c
        ../src/options.c
        ../src/win_compat.c
        ../src/defaults.h
        ../src/aperture.c
        ../src/geostragg.c
        ../src/script_command.c
        ../src/script_session.c
        ../src/script_file.c
        ../src/script_generic.h
        ../src/calibration.c
        ../src/prob_dist.c
        ../src/idf2jbs.c
        ../src/idfelementparsers.c
        ../src/idfparse.c
        ../src/nuclear_stopping.c
        ../src/fit_params.c
        ../src/stop.c
        ../src/des.c
        ../src/simulation_workspace.c
        ../src/sim_reaction.c
        ../src/sim_calc_params.c
        ../src/histogram.c
        ../src/gsl_inline.c
        ../src/scatint.c
        ../src/simulation2idf.c
        "$<$<BOOL:${JABS_PLUGINS}>:../src/plugin.c>"
       main.cpp
       mainwindow.cpp mainwindow.h mainwindow.ui
       fitdialog.cpp fitdialog.h fitdialog.ui
       qcustomplot.cpp qcustomplot.h
       spectrumplot.cpp spectrumplot.h
       codeeditor.cpp codeeditor.h
       highlighter.cpp highlighter.h
       plotdialog.h plotdialog.cpp plotdialog.ui
       preferencesdialog.h preferencesdialog.cpp preferencesdialog.ui
       icons.qrc
)

if (WIN32)
   set(app_icon_resource_windows "${CMAKE_CURRENT_SOURCE_DIR}/icon.rc")
   qt_add_executable(qjabs
	WIN32
        ${PROJECT_SOURCES}
	${app_icon_resource_windows}
    )
elseif(APPLE)
    set(MACOSX_BUNDLE_ICON_FILE icon.icns)
    set(app_icon_macos "${CMAKE_CURRENT_SOURCE_DIR}/icons/icon.icns")
    set_source_files_properties(${app_icon_macos} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources")
    qt_add_executable(qjabs MACOSX_BUNDLE
    ${PROJECT_SOURCES}
    ${app_icon_macos}
    )
    set_target_properties(qjabs PROPERTIES
        MACOSX_BUNDLE_BUNDLE_NAME JaBS
        MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
        MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
        MACOSX_BUNDLE_LONG_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}
        MACOSX_BUNDLE_ICON_FILE icon.icns
    )
else()
    qt_add_executable(qjabs
     ${PROJECT_SOURCES}
     )
endif()


target_include_directories(qjabs
        PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/.."
        ${GETOPT_INCLUDE_DIR}
        )

target_link_libraries(qjabs
    PRIVATE Qt6::Widgets
    PRIVATE Qt6::PrintSupport
    PRIVATE Qt6::Svg
    PRIVATE jibal
    PRIVATE GSL::gsl
    PRIVATE LibXml2::LibXml2
    PRIVATE "$<$<BOOL:${UNIX}>:m>"
    ${GETOPT_LIBRARY}
    gitwatcher
    )

if(OpenMP_C_FOUND)
    target_link_libraries(qjabs PUBLIC OpenMP::OpenMP_C)
    if (MSVC)
        add_compile_options(/openmp)
    endif ()
endif()

install(TARGETS qjabs
    BUNDLE DESTINATION .
    RUNTIME DESTINATION bin)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    install(FILES jabs.desktop DESTINATION share/applications)
    install(FILES jabs.xml DESTINATION share/mime/packages) 
    install(FILES icons/jabs.svg DESTINATION share/icons/hicolor/scalable/apps)
endif()
