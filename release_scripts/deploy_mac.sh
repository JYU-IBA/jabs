#!/bin/bash
cd .. #Everything relative to one directory up
#This script assumes you have installed JIBAL in ("install" directory)
install_dir="$(pwd)/install/"
jibal_data_file="jibal_data.tar.gz"
jibal_data_dir="$(pwd)/jibal_data/"
echo "Assuming JIBAL and QJaBS have been installed to $install_dir."
if [ ! -d "${install_dir}" ]; then
    echo "Install dir $install_dir} does not exist."
    exit 1;
fi

export PATH="$PATH:${install_dir}/bin"

if [ -f "${jibal_data_file}" ]; then
    echo "JIBAL data file exists."
else
    echo "Downloading JIBAL data file."
    curl "http://users.jyu.fi/~jaakjuli/jibal/data/data.tar.gz" -o "jibal_data.tar.gz"
fi
if echo "799a74a37cfd3d003a8b4914e6ff2c554c1b78eacbba84e3a6a414f982ce0bce  jibal_data.tar.gz" | shasum -a 256 --check --status; then
        echo "$jibal_data_file passes SHA sum check."
    else
        echo "$jibal_data_file does not pass SHA sum check."
        shasum -a 256 "$jibal_data_file"
        exit 1;
fi

rm -rf jibal_data
mkdir -p jibal_data
tar zxvf jibal_data.tar.gz -C "${jibal_data_dir}"

read version < version.txt
arch=$(uname -m)

echo "datadir = .
masses_file = masses.dat
abundances_file = abundances.dat
files_file = files.txt
assignments_file = assignments.txt
" > "${install_dir}/qjabs.app/Contents/Resources/jibal.conf"
cp "${install_dir}/share/jibal/masses.dat" "${install_dir}/share/jibal/abundances.dat" "${install_dir}/qjabs.app/Contents/Resources/"

echo "srim2013,srim2013.ele
yang,yang.stg
chu,chu.stg
bohr,bohr.stg
" > "${install_dir}/qjabs.app/Contents/Resources/files.txt"
cp "$jibal_data_dir"/*.ele "${install_dir}/qjabs.app/Contents/Resources/"
cp "$jibal_data_dir"/*.stg "${install_dir}/qjabs.app/Contents/Resources/"
cd "$install_dir"
mv qjabs.app JaBS.app
macdeployqt JaBS.app -dmg
mv JaBS.dmg "JaBS $version macOS $arch.dmg"
