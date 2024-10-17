REM This is mostly for deploying a Windows build on GitHub
REM Directories relative to one directory up (root of jabs repository)
REM This will assume you have run CMake install with --prefix=..\install (for JIBAL, JaBS and QJaBS) AND that qjabs has been built in ..\build_qjabs\ build directory (as a Release build)
cd ..
set INSTALL_DIR=..\install
set BUILD_DIR=..\build_qjabs\Release

curl "http://users.jyu.fi/~jaakjuli/jibal/data/data.tar.gz" -o jibal_data.tar.gz
mkdir "%BUILD_DIR%\JIBAL_data"
tar -zxvf jibal_data.tar.gz -C "%BUILD_DIR%\JIBAL_data"
copy "%INSTALL_DIR%\bin\*.dll" "%BUILD_DIR%"
copy "%INSTALL_DIR%\bin\*.exe" "%BUILD_DIR%"
del "%BUILD_DIR%\jibal.conf"
( 
  echo datadir = JIBAL_data
  echo masses_file = JIBAL_data\masses.dat
  echo abundances_file = JIBAL_data\abundances.dat
  echo files_file = JIBAL_data\files.txt
  echo assignments_file = JIBAL_data\assignments.txt
) > "%BUILD_DIR%\jibal.conf"
mkdir "%BUILD_DIR%\JIBAL_data"
copy "%INSTALL_DIR%\share\jibal\*.dat" "%BUILD_DIR%\JIBAL_data"
mkdir "%BUILD_DIR%\example"
copy example\experimental.dat "%BUILD_DIR%\example"
copy example\sample.txt "%BUILD_DIR%\example"
copy example\example.jbs "%BUILD_DIR%\example"
copy example\detector.jbs "%BUILD_DIR%\example"
mkdir "%BUILD_DIR%\example\tests"
copy example\tests\*.jbs "%BUILD_DIR%\example\tests"
copy example\tests\*.dat "%BUILD_DIR%\example\tests"
copy example\tests\*_ref.csv "%BUILD_DIR%\example\tests"
copy example\tests\*.r33 "%BUILD_DIR%\example\tests"
copy example\tests\*.txt "%BUILD_DIR%\example\tests"
copy LICENSE.txt "%BUILD_DIR%"
copy README.md "%BUILD_DIR%"
copy CITATION.cff "%BUILD_DIR%"
copy version.txt "%BUILD_DIR%"
(
  echo Windows release compiled %date% %time% by %COMPUTERNAME%
) > "%BUILD_DIR%/release.txt"
ver >> "%BUILD_DIR%/release.txt"
git describe --tags --dirty --broken --long --always >> "%BUILD_DIR%/release.txt"
