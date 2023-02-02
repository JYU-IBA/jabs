#!/bin/bash

jabslogfile="tests.log"

./clean.sh

error_exit() {
    echo "Error when running test ${1}! Aborted. See log in $jabslogfile.";
    exit 1;
}

for jabsfile in *.jbs; do
    echo "Running JaBS script $jabsfile";
#    if [ "$jabsfile" == "ds.jbs" ]; then    
#        echo "Skipping DS test."; # Takes too long...
#        continue
#    fi
#    if [ "$jabsfile" == "cornell_complicated.jbs" ]; then
#        echo "Skipping cornell_complicated test."; # Has DS too
#        continue
#    fi
    if [ "$jabsfile" == "plugin_cs.jbs" ]; then
        echo "Skipping plugin test."; # JaBS can be compiled without plugin support 
        continue
    fi      
    jabs "$jabsfile" 2>>"$jabslogfile" ||  error_exit "$jabsfile"
done

echo "All tests passed! Output saved in $jabslogfile."

#TODO: inspect generated CSV files
