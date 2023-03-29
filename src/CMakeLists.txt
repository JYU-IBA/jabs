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

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_C_STANDARD 99)

configure_file(version.h.in version.h @ONLY)
find_package(GSL 2.6 REQUIRED)
find_package(Jibal 0.3.10 REQUIRED)
find_package(LibXml2 REQUIRED)
find_package(OpenMP)

if(WIN32)
    find_path(GETOPT_INCLUDE_DIR getopt.h)
    find_library(GETOPT_LIBRARY getopt)
endif()

add_executable(jabs
        main.c
        sample.c brick.c ion.c simulation.c reaction.c options.c
        spectrum.c fit.c fit_params.c rotate.c detector.c jabs.c
        roughness.c  script.c  generic.c message.c aperture.c
        geostragg.c  script_command.c script_session.c script_file.c
        calibration.c prob_dist.c idf2jbs.c idfelementparsers.c
        idfparse.c nuclear_stopping.c stop.c des.c
        simulation_workspace.c sim_reaction.c sim_calc_params.c
        histogram.c gsl_inline.c scatint.c
        "$<$<BOOL:${JABS_PLUGINS}>:plugin.c>"
        "$<$<BOOL:${WIN32}>:win_compat.c>"
        )

add_executable(amsel amsel.c)
target_link_libraries(amsel PRIVATE jibal)

target_include_directories(jabs PRIVATE ${GETOPT_INCLUDE_DIR})

target_link_libraries(jabs
    PRIVATE jibal
    PRIVATE GSL::gsl
    PRIVATE LibXml2::LibXml2
    PRIVATE "$<$<BOOL:${UNIX}>:m>"
    ${GETOPT_LIBRARY} gitwatcher
    )


if(OpenMP_C_FOUND)
    target_link_libraries(jabs PUBLIC OpenMP::OpenMP_C)
endif()


install(TARGETS jabs RUNTIME DESTINATION bin)

if(0)
    if(DEBUG_MODE)
        add_subdirectory(des_table_test)
    endif()
endif()
