# Binary distribution

Standalone Windows and macOS disk images (dmg) are found in [GitHub releases](https://github.com/JYU-IBA/jabs/releases).

# Compiling JaBS from sources

JaBS uses a standard CMake build system. See instructions on [running CMake](https://cmake.org/runningcmake/) or follow instructions below.

## Microsoft Windows 10:

1. Install [JIBAL](https://github.com/JYU-IBA/jibal/blob/master/INSTALL.md).
2. Install *libxml2* using vcpkg (follow the JIBAL instructions mutatis mutandis)
3. See and run [deploy_win.bat](release_scripts/deploy_win.bat) file to prepare a standalone installation. There are hard coded locations in the file, so be sure to change those if necessary.
## Linux
1. Install [JIBAL](https://github.com/JYU-IBA/jibal/blob/master/INSTALL.md).
2. Install *libxml2* using your distributions package manager
    - On Ubuntu / Debian / Raspberry Pi OS: `apt install libxml2-dev`
    - On Arch: `pacman -S libxml2`
2. Run the following:

        $ git clone https://github.com/JYU-IBA/jabs.git
        $ mkdir build && cd build
        $ cmake ../
        $ make
        $ sudo make install

3. Compiling Qt GUI, install Qt 6 and try similar build steps in [qjabs](qjabs/) directory

## MacOS:
1. Install [Homebrew](https://brew.sh/)
2. Run the following to install command line version of JaBS with all dependencies:
    
        $ brew tap JYU-IBA/iba
        $ brew install jabs-cli

4. If you want to develop JaBS and not just use it, follow Linux instructions above, install *libxml2* and *qt6* using Homebrew.
5. Disk image (dmg) can be created using [deploy_mac.sh](release_scripts/deploy_mac.sh) script, assuming JIBAL data can be found from home directory.
