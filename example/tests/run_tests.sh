#!/bin/bash
for jabsfile in *.jbs; do
    echo "Running JaBS script $jabsfile";
    if [ "$jabsfile" == "ds.jbs" ]; then
        echo "Skipping test."; # Takes too long...
    else
        jabs "$jabsfile"|| exit 1;
    fi
done

#TODO: inspect generated CSV files
