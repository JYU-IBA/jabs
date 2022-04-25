#!/bin/bash
build_dir="build-qjabs-Desktop_x86_darwin_generic_mach_o_64bit-Release"
jibal_install_prefix="/usr/local"
read version < version.txt
arch=$(uname -m)
rm -f example.zip
cd "$build_dir" || exit 1
make
echo "datadir = .
masses_file = masses.dat
abundances_file = abundances.dat
files_file = files.txt
assignments_file = assignments.txt
" > "qjabs.app/Contents/Resources/jibal.conf"
cp "$jibal_install_prefix/share/jibal/masses.dat" "$jibal_install_prefix/share/jibal/abundances.dat" "qjabs.app/Contents/Resources/"

echo "srim2013,srim2013.ele
yang,yang.stg
chu,chu.stg
bohr,bohr.stg
" > "qjabs.app/Contents/Resources/files.txt"
cp "$HOME/.jibal/srim2013.ele" "$HOME/.jibal/yang.stg" "$HOME/.jibal/chu.stg" "$HOME/.jibal/bohr.stg" "qjabs.app/Contents/Resources/"
mv qjabs.app JaBS.app
macdeployqt JaBS.app -dmg
rm -rf JaBS.app
mv JaBS.dmg "../JaBS $version macOS $arch.dmg"
cd ..
find example|grep -v -e "/\." -e "_out\."|zip example -@
