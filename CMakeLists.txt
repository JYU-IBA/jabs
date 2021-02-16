cmake_minimum_required(VERSION 3.15)
project(jabs
        VERSION 0.1.0
        DESCRIPTION "Jaakko's Backscattering Simulator (JaBS)"
        LANGUAGES C
        )

set(CMAKE_C_STANDARD 99)

configure_file(defaults.h.in defaults.h @ONLY)
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

add_executable(jabs main.c sample.c sample.h brick.c ion.c ion.h simulation.c simulation.h reaction.c reaction.h)
target_include_directories(jabs PRIVATE ${GETOPT_INCLUDE_DIR})

target_link_libraries(jabs
    PRIVATE Jibal 
    PRIVATE GSL::gsl
    PRIVATE "$<$<BOOL:${UNIX}>:m>" ${GETOPT_LIBRARY}
    )

install(TARGETS jabs RUNTIME DESTINATION bin)
