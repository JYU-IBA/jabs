#!/bin/bash
cd .. #Everything relative to one directory up
#This script assumes you have installed JIBAL somewhere (so that CMake can find it and jibaltool can be run)
#Try running build_jibal.sh and build_qjabs.sh first
builddir="../build_qjabs"
resources="${builddir}/qjabs.app/Contents/Resources/"
if ! jibaltool status; then echo
    "Install JIBAL first. CMake should be able to find it and I need to run jibaltool."
    "Check or run build_jibal.sh"
fi
jibal_datadir="$(jibaltool config 2>&1|grep "^datadir = "|sed "s/^datadir = //")"
jibal_gsto_data_file="jibal_data.tar.gz"
jibal_gsto_data_dir="$(pwd)/jibal_data/"
echo "Assuming JIBAL data (masses.dat etc) is found in $jibal_datadir."

echo "Assuming QJaBS is built in $builddir"
if [ ! -d "${builddir}" ]; then
    echo "Build directory ${builddir} does not exist."
    echo "Try running build_qjabs.sh"
    exit 1;
fi


read version < version.txt
arch=$(uname -m)
echo "Arch $arch"
output="JaBS $version macOS $arch.dmg"

if [ -f "${jibal_gsto_data_file}" ]; then
    echo "JIBAL data file exists."
else
    echo "Downloading JIBAL data file."
    curl "http://users.jyu.fi/~jaakjuli/jibal/data/data.tar.gz" -o "jibal_data.tar.gz"
fi
if echo "799a74a37cfd3d003a8b4914e6ff2c554c1b78eacbba84e3a6a414f982ce0bce  jibal_data.tar.gz" | shasum -a 256 --check --status; then
        echo "$jibal_gsto_data_file passes SHA sum check."
    else
        echo "$jibal_gsto_data_file does not pass SHA sum check."
        shasum -a 256 "$jibal_gsto_data_file"
        exit 1;
fi

rm -rf jibal_data
rm -f "$output"
mkdir -p jibal_data
tar zxvf jibal_data.tar.gz -C "${jibal_gsto_data_dir}"

echo "datadir = .
masses_file = masses.dat
abundances_file = abundances.dat
files_file = files.txt
assignments_file = assignments.txt
" > "${resources}/jibal.conf"
cp "${jibal_datadir}/masses.dat" "${jibal_datadir}/abundances.dat" "${resources}"
echo "srim2013,srim2013.ele
yang,yang.stg
chu,chu.stg
bohr,bohr.stg
" > "${resources}/files.txt"
cp "$jibal_gsto_data_dir"/*.ele "${resources}"
cp "$jibal_gsto_data_dir"/*.stg "${resources}"
cp "LICENSE.txt" "${resources}"
cd "$builddir"
echo "Operating in $pwd"
#mv qjabs.app JaBS.app
rm -rf JaBS.app
cp -R qjabs.app JaBS.app
echo "Running macdeployqt"
if [ -z "${JABS_SIGN_IDENTITY}" ]; then 
    codesign_opt="";
else
    echo "Signing with $JABS_SIGN_IDENTITY"
    codesign_opt=-codesign="${JABS_SIGN_IDENTITY}";
fi

#cp lib/libjibal.0.dylib JaBS.app/Contents/Frameworks/
#install_name_tool -change /opt/homebrew/opt/gsl/lib/ "@executable_path/../Frameworks/" JaBS.app/Contents/Frameworks/libjibal.0.dylib
macdeployqt6 JaBS.app -dmg -no-strip ${codesign_opt}
mv JaBS.dmg "$output"
