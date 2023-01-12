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
