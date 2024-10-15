#!/bin/bash
cd .. #Everything relative to root of repository
jibal_src_dir="../jibal"
builddir="../build_jibal"
if [ ! -d "$jibal_src_dir" ]; then
    echo "Tried to find JIBAL source code from directory $jibal_src_dir (relative to $(pwd))"
    echo "Perhaps run this in the parent directory (of $(pwd)): "
    echo "git clone https://github.com/JYU-IBA/jibal.git"
    exit 1
fi
mkdir -p "$builddir";
if ! cmake -S "$jibal_src_dir" -B "$builddir" -DCMAKE_BUILD_TYPE=Release; then
    echo "Could not configure using CMake"
    exit 1;
fi

if ! cmake --build "$builddir" --config Release; then
    echo "Could not build."
    exit 1;
fi

if ! sudo cmake --install "$builddir"; then
    echo "Could not install."
    exit 1;
fi
