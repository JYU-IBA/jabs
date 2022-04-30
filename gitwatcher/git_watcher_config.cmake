set(PRE_CONFIGURE_FILE "${CMAKE_CURRENT_LIST_DIR}/git.c.in")
set(POST_CONFIGURE_FILE "${CMAKE_CURRENT_BINARY_DIR}/git.c")
set(GIT_FAIL_IF_NONZERO_EXIT FALSE)
include(${CMAKE_CURRENT_LIST_DIR}/git_watcher.cmake)

# Create a library out of the compiled post-configure file.
add_library(gitwatcher STATIC ${POST_CONFIGURE_FILE})
target_include_directories(gitwatcher PUBLIC ${CMAKE_CURRENT_LIST_DIR})
add_dependencies(gitwatcher check_git)
