cmake_minimum_required(VERSION 3.15)
project(idf2jbs
    VERSION 0.0.1
    DESCRIPTION "IDF to JaBS converter"
    LANGUAGES C
)

option(DEBUG_MODE "Debug mode" OFF)
if(DEBUG_MODE)
    if (MSVC)
        add_compile_options(/W4 /WX /DDEBUG)
    else()
        add_compile_options(-g -O0 -Wall -Wextra -pedantic -DDEBUG)
    endif()
endif()

find_package(LibXml2 REQUIRED)
find_package(Jibal REQUIRED)

add_executable(idf2jbs main.c idf2jbs.c)

target_link_libraries(idf2jbs PRIVATE LibXml2::LibXml2 PRIVATE jibal)

install(TARGETS idf2jbs RUNTIME DESTINATION bin)
