cmake_minimum_required(VERSION 3.15)
project(simu
        VERSION 0.0.1
        DESCRIPTION "RBS simulator. Shitty."
        LANGUAGES C
        )

find_package(Jibal REQUIRED)
option(DEBUG_OUTPUT_ENABLE "Enable debug output" OFF)
if(DEBUG_OUTPUT_ENABLE)
    add_compile_options(-DDEBUG)
endif()
add_executable(simu main.c)
target_link_libraries(simu Jibal)

install(TARGETS simu RUNTIME DESTINATION bin)
