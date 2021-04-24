#!/usr/bin/env gnuplot
set title "JABS Example, 2 MeV 4He, thin Au on SiO2/Si"
set datafile separator ','
set key autotitle columnhead
file='out.csv'

set xrange [200:2000]

set xlabel "Energy (keV)"
set ylabel "Counts"

plot \
file u 4:3 w histeps lc rgbcolor "blue", \
file u 4:2 w lines lc rgbcolor "black",\
#file u 1:2 w steps lc rgbcolor "black",\
