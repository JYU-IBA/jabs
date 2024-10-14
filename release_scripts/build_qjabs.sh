#!/bin/bash
cd ..
builddir="build_qjabs"
mkdir -p "$builddir"
if ! cmake --fresh -S qjabs -B "$builddir" -DCMAKE_BUILD_TYPE=Release; then
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
