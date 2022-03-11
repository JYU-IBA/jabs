cmake_minimum_required(VERSION 3.19)
if(WIN32 AND DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
            CACHE STRING "")
    message(STATUS "CMake toolchain file ${CMAKE_TOOLCHAIN_FILE}")
endif()
project(qjabs
    VERSION 0.5.3
    DESCRIPTION "Jaakko's Backscattering Simulator (JaBS) for Qt"
    LANGUAGES C CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MACOSX_BUNDLE TRUE)

configure_file(version.h.in version.h @ONLY)

find_package(Qt6 COMPONENTS Widgets REQUIRED PrintSupport REQUIRED)

if(${QT_VERSION_MAJOR} NOT EQUAL 6)
    message("Qt major version is ${QT_VERSION_MAJOR}. This program is only for Qt6.")
endif()

find_package(GSL 2.6 REQUIRED)
find_package(Jibal 0.3.5 REQUIRED)
if(WIN32)
    find_path(GETOPT_INCLUDE_DIR getopt.h)
    find_library(GETOPT_LIBRARY getopt)
endif()

set(PROJECT_SOURCES
        ../sample.c
        ../sample.h
        ../brick.c
        ../ion.c
        ../ion.h
        ../simulation.c
        ../simulation.h
        ../reaction.c
        ../reaction.h
       ../sample.c
       ../sample.h
       ../spectrum.c
       ../spectrum.h
       ../fit.c
       ../fit.h
       ../rotate.c
       ../rotate.h
       ../jabs.c
       ../jabs.h
       ../generic.c
       ../generic.h
       ../detector.c
       ../detector.h
       ../roughness.c
       ../roughness.h
       ../script.c
       ../script.h
       ../jabs.c
       ../jabs.h
       ../options.c
       ../options.h
       ../message.h
        ../win_compat.c
        ../win_compat.h
        ../defaults.h
        ../aperture.c
        ../aperture.h
        ../geostragg.c
        ../geostragg.h
        ../script_command.c
        ../script_command.h
        ../script_session.c
        ../script_session.h
        ../script_file.c
        ../script_file.h
        ../script_generic.h
        ../calibration.c
        ../calibration.h
       main.cpp
       mainwindow.cpp
       mainwindow.h
       mainwindow.ui
       qcustomplot.cpp
       qcustomplot.h
       spectrumplot.cpp
       spectrumplot.h
       icons.qrc
       codeeditor.cpp
       codeeditor.h
       highlighter.cpp
       highlighter.h
)

set(MACOSX_BUNDLE_ICON_FILE icon.icns)
set(app_icon_macos "${CMAKE_CURRENT_SOURCE_DIR}/icons/icon.icns")
set_source_files_properties(${app_icon_macos} PROPERTIES
       MACOSX_PACKAGE_LOCATION "Resources")

if (WIN32)
#   set(app_icon_resource_windows "${CMAKE_CURRENT_SOURCE_DIR}/icon.rc")
#    qt_add_executable(qjabs ${PROJECT_SOURCES} ${app_icon_resource_windows})
    qt_add_executable(qjabs
        ${PROJECT_SOURCES}
    )
elseif(APPLE)
   qt_add_executable(qjabs
    MACOSX_BUNDLE
    ${PROJECT_SOURCES}
    ${app_icon_macos}
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
    PRIVATE jibal
    PRIVATE GSL::gsl
    PRIVATE "$<$<BOOL:${UNIX}>:m>"
    ${GETOPT_LIBRARY}
    )

set_target_properties(qjabs PROPERTIES
    MACOSX_BUNDLE_BUNDLE_NAME QJaBS
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE_ICON_FILE icon.icns
)
