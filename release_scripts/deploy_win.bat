REM This script will make a standalone JaBS installation in qjabs\build\Release
REM Install JIBAL first and run the bootstrap tool
set JIBAL_DIR=C:\Program Files\jibal
REM Qt6 must be installed too
set QT6_DIR=C:\Qt\6.4.2
REM And vcpkg
set VCPKG_DIR=C:\vcpkg
REM Try to remove CMakeCache.txt files if your configuration (e.g. compiler) changes or you want to do a fresh build
cd ..
mkdir build
cd build
REM del CMakeCache.txt
cmake -G "Visual Studio 16 2019" -A x64 -DCMAKE_TOOLCHAIN_FILE="%VCPKG_DIR%\scripts\buildsystems\vcpkg.cmake" ../
cmake --build . --target ALL_BUILD --config Release
cd ..

cd qjabs
mkdir build
cd build
REM del CMakeCache.txt
cmake -DCMAKE_PREFIX_PATH="%QT6_DIR%\msvc2019_64\lib\cmake" -G "Visual Studio 16 2019" -A x64 -DCMAKE_TOOLCHAIN_FILE="%VCPKG_DIR%\scripts\buildsystems\vcpkg.cmake" ../
cmake --build . --target ALL_BUILD --config Release
REM Windeployqt will handle most dlls and other Qt dependencies
%QT6_DIR%\msvc2019_64\bin\windeployqt.exe Release\qjabs.exe
copy ..\..\build\src\Release\jabs.exe Release
copy "%JIBAL_DIR%\bin\*.dll" Release
del Release\jibal.conf
( 
  echo datadir = JIBAL_data
  echo masses_file = JIBAL_data\masses.dat
  echo abundances_file = JIBAL_data\abundances.dat
  echo files_file = JIBAL_data\files.txt
  echo assignments_file = JIBAL_data\assignments.txt
) > Release\jibal.conf
mkdir Release\JIBAL_data
copy "%JIBAL_DIR%\share\jibal\*.dat" Release\JIBAL_data
copy "%AppData%\Jibal\*.ele" Release\JIBAL_data
copy "%AppData%\Jibal\*.stg" Release\JIBAL_data
copy "%AppData%\Jibal\files.txt" Release\JIBAL_data
copy "%AppData%\Jibal\assignments.txt" Release\JIBAL_data
mkdir Release\example
copy "..\..\example\experimental.dat" Release\example
copy "..\..\example\sample.txt" Release\example
copy "..\..\example\example.jbs" Release\example
copy "..\..\example\detector.jbs" Release\example
mkdir Release\example\tests
copy ..\..\example\tests\*.jbs Release\example\tests\
copy ..\..\example\tests\*.dat Release\example\tests\
copy ..\..\example\tests\*_ref.csv Release\example\tests\
copy ..\..\example\tests\*.r33 Release\example\tests\
copy ..\..\example\tests\*.txt Release\example\tests\
cd ..
cd ..
