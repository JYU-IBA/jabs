#!/bin/bash
#Compiles the command line interface version of JaBS
cd .. #Every path relative to root of repository
builddir="../build_jabs"
mkdir -p "$builddir"
if ! cmake --fresh -S . -B "$builddir" -DCMAKE_BUILD_TYPE=Release; then
    echo "Could not configure using CMake"
    exit 1;
fi

if ! cmake --build "$builddir" --config Release; then
    echo "Could not build."
    exit 1;
fi

#if ! cmake --install "$builddir"; then
#    echo "Could not install."
#    exit 1;
#fi
