#One day after "make install" everything worked on macOS. The next day dyld complains "Symbol not found" when trying to run an executable that uses JIBAL and is installed in /usr/local/bin. This made people angry and sad.
# 1. JIBAL dylib does have the symbols. 
# 2. Old (previously compiled) executables work, new ones don't. Same JIBAL library is installed, old versions have been removed.
# 3. "otool -L" says "@rpath/libjibal.0.dylib" for both old and new executables.
# 4. "otool -l" doesn't list any LC_RPATHs for the old executables (that still work) nor for the new ones (that don't work).
# 5. CMake is updated to latest version (3.27.5). No difference. Apple clang version 15.0.0 (clang-1500.0.40.1) Target: x86_64-apple-darwin22.6.0
#
# So after a lot of googling we have a temporary fix. This file. Hopefully this fix doesn't break anything else. It should essentially just put the installation prefix (e.g. /usr/local) to rpath in some cases. Like in my case.

#The following is directly from https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling (accessed Sep 19 2023)

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")
