#!/bin/bash
if builtin type -P "jabs" >/dev/null; then
    clear;
    echo -e "I am going to run a JaBS fit with these arguments:\n\n---arguments.txt---\n"
    xargs echo < "arguments.txt"
    echo -e "\n-------------------\n\n"
    echo "See out.csv for simulated and experimental data.";
    echo "Sleeping for 10 seconds to give you time to read this."
    echo "Note that the fitting will take a while, there is roughness..."
    sleep 10;
    xargs jabs < "arguments.txt" 
    sleep 2;
    echo ""
    echo "Thanks for running the example script! You can use it as your template."
    echo "Just modify arguments.txt to suit your needs."
    echo "See output in out.csv, there is a Gnuplot script to plot it in plot.plt."
    echo "Fitted sample and detector models are saved in sample_fitted.txt and detector_fitted.txt."
    echo "These files could be used as input for another fit."
else
    echo "Install jabs first. See instructions in README.md";
fi
