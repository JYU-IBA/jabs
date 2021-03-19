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

        $ mkdir build
        $ cd build
        $ cmake ../
        $ make
       
## Usage

Several parameters can be set from the command line. See `jabs -h` and try something like this:

~~~~
$ ./jabs -E 2MeV --alpha=10deg --theta=170deg --out=spectrum.csv Au 500tfu SiO2 1000tfu Si 10000tfu
~~~~


## Fitting

There is an experimental fitting feature too, activated with the `--fit` or `-F` option. The multidimensional nonlinear least-squares fitting is based on [GSL multifit](https://www.gnu.org/software/gsl/doc/html/nls.html).


The user can currently give only a single range to fit. This is specified with `--fit_low` and `--fit_high`. Default values are 10% of highest channel and the highest channel, respectively.

Use `--fit_vars` to give a comma separate list of variables to fit. Currently supported variables are `slope,offset,resolution` related to detector parameters, 
`calib`which is equivalent to the previous triplet, `fluence`, and `thicknessN` where `N` is number of a layer, numbering starts from 1 (surface).

Example: `jabs --fit_vars=calib,fluence,thickness1 --fit_low=250 ...`
