#!/usr/bin/env gnuplot
#This is a Gnuplot example of how to plot JaBS CSV output
set title "JABS Example, 2 MeV 4He, thin Au on SiO2/Si"
set datafile separator ','
set key autotitle columnhead
file='example_out.csv'

set xrange [200:2000]

set xlabel "Energy (keV)"
set ylabel "Counts"

plot \
file u 2:4 w histeps lc rgbcolor "blue",\
file u 2:3 w lines lc rgbcolor "black",\
