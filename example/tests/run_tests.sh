#!/bin/bash

error_exit() {
    echo "Error when running test ${1}! Aborted.";
    exit 1;
}
for jabsfile in *.jbs; do
    echo "Running JaBS script $jabsfile";
    if [ "$jabsfile" == "ds.jbs" ]; then
        echo "Skipping DS test."; # Takes too long...
    elif [ "$jabsfile" == "plugin_cs.jbs" ]; then
        echo "Skipping plugin test."; # JaBS can be compiled without plugin support 
    else
        jabs "$jabsfile"||  error_exit "$jabsfile"
    fi
done

#TODO: inspect generated CSV files
