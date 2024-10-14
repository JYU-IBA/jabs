# Binary distribution

Standalone Windows and macOS disk images (dmg) are found in [GitHub releases](https://github.com/JYU-IBA/jabs/releases).

# Compiling JaBS from sources

JaBS uses a standard CMake build system. See instructions on [running CMake](https://cmake.org/runningcmake/) or follow instructions below.

JaBS depends on [GSL](https://www.gnu.org/software/gsl/), [JIBAL](https://github.com/JYU-IBA/jibal), [libxml2](https://gitlab.gnome.org/GNOME/libxml2/-/wikis/home), and the GUI (QJaBS)
 depends additionally on [Qt 6](https://www.qt.io/) and [QCustomPlot](https://www.qcustomplot.com/).

Preferred compiler on Windows is MSVC 2019 (for Qt compatibility), Clang on macOS and GCC on Linux.

Creating binary distributions may be challenging, see [GitHub workflows](.github/workflows). The Windows binary distribution is created with this workflow.

There are scripts that manage downloading QCustomPlot and JIBAL datafiles (for Windows and macOS).
## Linux
1. Install dependencies using your distributions package manager
   - On Ubuntu / Debian / Raspberry Pi OS: `sudo apt install  libgsl27 libgsl-dev libxml2 libxml2-dev libreadline-dev`
   - On Arch: `sudo pacman -S gsl libxml2 readline`
2. Obtain source codes of [JIBAL](https://github.com/JYU-IBA/jibal/) and JaBS (this repository)
   ```shell
   git clone https://github.com/JYU-IBA/jibal.git
   git clone https://github.com/JYU-IBA/jabs.git
   ```
3. Compile JIBAL and install it to /usr/local (recommended, but you may also install it to some other prefix)
   ```shell
   cmake -B jibal/build_jibal -S jibal
   cmake --build jibal/build_jibal
   sudo cmake --install build_jibal
   ```

4. Compile and install JaBS (CLI version)
   ```shell
   cmake -B jabs/build_jabs -S jabs
   cmake --build jabs/build_jabs
   sudo cmake --install jabs/build_jabs
   ```
5. Install dependencies for the Qt version:
   ```shell
   sudo apt install qt6-base-dev
   ```
6. Run script `get_external_files.sh` or download [QCustomPlot](https://www.qcustomplot.com/index.php/download) source code manually. The cpp and h files must be placed in qjabs directory.

7. Compile and install QJaBS (Qt GUI version)
   ```shell
   cmake -B jabs/build_qjabs -S jabs/qjabs
   cmake --build jabs/build_qjabs
   sudo cmake --install jabs/build_qjabs
   ```

If you get errors when running JaBS, check that the JIBAL library can be found, e.g.
```shell
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib"
sudo ldconfig
```

There are also build scripts in the [release scripts](release_scripts/) directory. See what they do before running them.

You will also need some stopping data for JIBAL:
```shell
curl "http://users.jyu.fi/~jaakjuli/jibal/data/data.tar.gz" -o "jibal_data.tar.gz"
tar zxvf jibal_data.tar.gz -C "$HOME/.jibal"
```

## Microsoft Windows 10:

1. Install Build tools for [Visual Studio 2019.](https://visualstudio.microsoft.com/downloads/)
2. Install Qt 6.7.2 (or maybe something later). This might be easier to do with [unofficial tools](https://github.com/miurahr/aqtinstall/).
3. Install *gsl*, *getopt*, *libxml2* using vcpkg. Make sure your triplet is x64-windows-release.
4. Follow steps for Linux above, but
   - Run the first cmake (configure) with: `cmake -DCMAKE_TOOLCHAIN_FILE=C:/vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows-release -DVCPKG_MANIFEST_MODE=OFF -DCMAKE_BUILD_TYPE=Release -B build_dir -S source_dir`
   - Install JIBAL (and everything else) to some prefix (here: *install*) with `cmake --install build_jibal --prefix install`
   - Give subsequent JaBS and QJaBS CMake configure steps `-DCMAKE_PREFIX_PATH=install` so that JaBS can find the installed JIBAL library from the directory given by the prefix
5. Run [Qt Windows Deployment tool](https://doc.qt.io/qt-6/windows-deployment.html) (windeployqt) on qjabs.exe:
   ```shell
       windeployqt build_qjabs\Release\qjabs.exe
   ```
6. See and run [deploy_win.bat](release_scripts/deploy_win.bat) file to prepare a standalone installation. This is mostly intended to be run by GitHub, so some hard-coded directories may be different than in the instructions above.

## MacOS:
1. Install [Homebrew](https://brew.sh/)
2. Run the following to install command line version of JaBS with all dependencies:
    
        $ brew tap JYU-IBA/iba
        $ brew install jabs-cli

3. If you want to develop JaBS and not just use it, follow Linux instructions above, install *libxml2*, *readline* and *qt6* using Homebrew. Installing *libomp* from homebrew is also recommended to enable OpenMP parallel processing. You may need to set the *OpenMP_ROOT* environment variable to get CMake to find it on Apple silicon, like this:

        $ export OpenMP_ROOT="/opt/homebrew/opt/libomp/"

4. Disk image (dmg) for QJaBS app can be created using [deploy_mac.sh](release_scripts/deploy_mac.sh) script. Note that some assumptions are made, so read it be fore using it. Running cmake --install is recommended for JIBAL and JaBS, but not for QJaBS.