REM This is mostly for deploying a Windows build on GitHub
cd ..
set JIBAL_DIR=.\install
set BUILD_DIR=.\build_qjabs\Release

curl "http://users.jyu.fi/~jaakjuli/jibal/data/data.tar.gz" -o jibal_data.tar.gz
mkdir "%BUILD_DIR%\JIBAL_data"
tar -zxvf jibal_data.tar.gz -C "%BUILD_DIR%\JIBAL_data"
copy "%JIBAL_DIR%\bin\*.dll" "%BUILD_DIR%"
del "%BUILD_DIR%\jibal.conf"
( 
  echo datadir = JIBAL_data
  echo masses_file = JIBAL_data\masses.dat
  echo abundances_file = JIBAL_data\abundances.dat
  echo files_file = JIBAL_data\files.txt
  echo assignments_file = JIBAL_data\assignments.txt
) > "%BUILD_DIR%\jibal.conf"
mkdir "%BUILD_DIR%\JIBAL_data"
copy "%JIBAL_DIR%\share\jibal\*.dat" Release\JIBAL_data
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
