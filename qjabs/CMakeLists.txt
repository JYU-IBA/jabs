cmake_minimum_required(VERSION 3.15)

project(qjabs
    VERSION 0.4.1
    DESCRIPTION "Jaakko's Backscattering Simulator (JaBS) for Qt"
    LANGUAGES C CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

configure_file(version.h.in version.h @ONLY)

find_package(Qt6 COMPONENTS Widgets REQUIRED PrintSupport REQUIRED)

if(${QT_VERSION_MAJOR} NOT EQUAL 6)
    message("Qt major version is ${QT_VERSION_MAJOR}. This program is only for Qt6.")
endif()

find_package(GSL 2.6 REQUIRED)
find_package(Jibal 0.3.3 REQUIRED)

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
       main.cpp
       mainwindow.cpp
       mainwindow.h
       mainwindow.ui
       qcustomplot.cpp
       qcustomplot.h
       spectrumplot.cpp
       spectrumplot.h
)

qt_add_executable(qjabs
    MANUAL_FINALIZATION
    ${PROJECT_SOURCES}
)

target_include_directories(qjabs PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")

target_link_libraries(qjabs
    PRIVATE Qt6::Widgets
    PRIVATE Qt6::PrintSupport
    PRIVATE jibal
    PRIVATE GSL::gsl
    PRIVATE "$<$<BOOL:${UNIX}>:m>"
    )

set_target_properties(qjabs PROPERTIES
    MACOSX_BUNDLE_GUI_IDENTIFIER jyu.fi
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
)

qt_finalize_executable(qjabs)
