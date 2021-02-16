# JaBS - Jaakko's Backscattering Simulator

Simulates RBS spectra rapidly.
    
    Copyright (C) 2021 Jaakko Julin
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

## Installation

1. Install [JIBAL](https://github.com/JYU-IBA/jibal). You need [GSL](https://www.gnu.org/software/gsl/) for both JIBAL and JaBS.
2. Clone this repository, e.g.

        $ git clone git@gitlab.jyu.fi:iba/jabs.git

3. Build, see [CMake instructions](https://cmake.org/runningcmake/), or try...

        mkdir build
        cd build
        cmake ../
        make
       
## Usage

Several (undocumented) parameters can be set from the command line. Try something like this:

~~~~
./jabs -f -I 4He -E 2MeV -o tmp.csv --alpha 10deg --beta 10deg --theta 160deg -R 20keV Au 500tfu SiO2 1000tfu Si 10000tfu
~~~~

