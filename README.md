# JaBS - Jaakko's Backscattering Simulator

Simulates RBS, ERD and NRA spectra rapidly.
    
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

Please see additional notes regarding copyrights of included open source code at the end of this document.

## Binary packages
Some ready-to-use [binary packages](http://users.jyu.fi/~jaakjuli/jabs/binaries/) for Windows and macOS may be available.

## Installation from sources

1. Install [JIBAL](https://github.com/JYU-IBA/jibal), preferably from sources. You'll need [GSL](https://www.gnu.org/software/gsl/) for both JIBAL and JaBS.
2. Clone this repository, e.g.

        $ git clone git@github.com:JYU-IBA/jabs.git


3. Build, see [CMake instructions](https://cmake.org/runningcmake/), or try...

        $ mkdir build
        $ cd build
        $ cmake ../
        $ make

## Graphical user interface

Simple Qt-based [GUI](qjabs/) is provided. Scripts can be edited and run, simulated spectra and experimental spectra are plotted (if available). Many GUI-specific features are missing, but the scripts themselves should work as described below. Currently some stopping forces and straggling data is not reloaded (from JIBAL) between runs and restarting the program may be necessary to reset these.

      
## Interactive or scripted usage

This is the preferred way of using JaBS. Some command line options (see below) can be used in conjunction with the interactive or scripted mode.

Launch JaBS in interactive mode simply by running `jabs` or `jabs --interactive`. Scripts (i.e. files with commands identical to interactive mode input) can be given as command line parameters e.g. `jabs script.jbs`, piped in `jabs < script.jbs` or run using the interactive mode. Filename from command line must not be `sample` (reserved keyword for giving the sample on the command line, see below).

The interactive mode should be self-explanatory and an internal help is provided. Please see [the example script](example/example.jbs) to get started.

## Command line usage

Several parameters can be set from the command line. See `jabs -h` and try something like this:

~~~~
$ ./jabs -E 2MeV --alpha=10deg --theta=170deg --out=spectrum.csv sample Au 500tfu SiO2 1000tfu Si 10000tfu
~~~~

Detector and sample can be read from files. The file formats are simple and human readable. Please see [the example](example).

## Features
### Implemented
 - Basic RBS spectrum simulation with atomic data, cross sections, electronic stopping and straggling given by [JIBAL](https://github.com/JYU-IBA/jibal)
 - Automatic recoil spectra when working in forward angles (ERDA)
 - Arbitrary cross sections and reactions from R33 files. Both EBS (Q-value is zero) as well as p-p NRA are implemented. 
 - Point-by-point and layered sample models
 - Roughness (gamma distribution)
 - Arbitrary geometry, detector and sample tilt can be expressed in arbitrary spherical coordinates
 - Linear detector calibration and constant detector energy resolution
 - Multilayer foil in front of detector
 - Reading of experimental data
 - Fitting to experimental spectrum (multidimensional non-linear least squares fitting)
 - Linear ad-hoc substrate channeling correction
 - Faster (less accurate) mode
 - User defined "molecules" i.e. elements with fixed concentration ratios (e.g. you can fit C in SiO2 without changing Si/O ratio)
 - Weighting of cross-sections by straggling
 - Dual scattering model (although it needs improvements before it is usable)
 - Simultaneous multi-detector simulation and fitting.
 - Stopping corrections can be supplied by user (Bragg correction) for a specific layer

### Not (yet) implemented, but planned
 - Support for more input and output data formats (CSV, IDF, ...)
 - More accurate handling of sharp peaks in cross sections (resonances). The current handling is quite accurate in most cases.
 - Multiple scattering (small angle)
 - Kinematic (geometric) broadening
 - Simulation of pile-up and dead time
 - Non-linear detector response and different response for different particles (e.g. alphas, protons)
 - Simulation of time-of-flight spectra
 - Publication quality plotting

### Distant future
 - Parallel processing of independent simulations
 - Advanced fitting algorithms
 - Fitting of spectra from different measurements (different beam, fluence etc for each simulation)

### Known issues
 - Occasional crashes when fitting, since some corner cases are not handled properly
 - Detector numbering and usability issues with multidetector mode 
 - Point-by-point profiles are not tested (but should work)
 - Ad-hoc channeling correction is the same for all detectors
 - Detector calibration is the same for all particles (issue with NRA and ERD, but not for RBS and EBS)
 - Dual scattering assumes first scattering is RBS (not ERD). Cross sections are not calculated accurately (must use integrated cross sections instead of approximating using differential cross sections since solid angles involved are large).
 - Dual scattering is benchmarked against SimNRA and is known to produce different results.
 - Adding detector related fit variables will add those for all detectors 

## Fitting

There is a fitting feature, activated with the `--fit` or `-F` option or `fit` script command. The multidimensional nonlinear least-squares fitting is based on [GSL multifit](https://www.gnu.org/software/gsl/doc/html/nls.html). The accuracy and sanity of fits must be evaluated by the user.

The user can give only a single range to fit from the command line. This is specified with `--fit_low` and `--fit_high`.

Use `--fit_vars` to give a comma separate list of variables to fit. Currently supported detector related variables are `slope,offset,resolution`
`calib` which is equivalent to the previous triplet and `solid`. Currently adding fit variables will add all them for each detector.

Beam `fluence` can be fitted. Don't try to combine this with `solid` for all detectors.

Layer thickness can be also fitted using variable `thicknessN` where `N` is number of a layer, numbering starts from 1 (surface). Roughness can be fitted similarly to thickness, use `roughN`. Concentrations can be fitted using `concN_I`, where `N` is the layer and `I` the element.

Example: `jabs --fit_vars=calib,fluence,thickness1 --fit_low=250 ...`

## Open source licenses

This repository contains following open source software, or parts of them:
 * Code derived from GSL (GNU Scientific Library) examples, e.g. in [fit.c](fit.c). Used under the terms of GPL v3.
 * QCustomPlot [qcustomplot.cpp](qjabs/qcustomplot.cpp) and [qcustomplot.h](qjabs/qcustomplot.h), Copyright (C) 2011-2021 Emanuel Eichhammer. Used under the terms of GPL v3.
 * Coordinate rotate routine from MCERD [rotate.c](rotate.c). Copyright Kai Arstila. User under the terms of GPL v2.
 * Some routines from NetBSD, Copyright (c) 2011 The NetBSD Foundation, Inc. See [win_compat.c](win_compat.c) for the full copyright notice.
 * Qt6 GUI is distributed under the provisions of LGPL version 3. Qt Toolkit is Copyright (C) 2018 The Qt Company Ltd. and other contributors.
 * Material Design [icons](qjabs/icons) licensed under [Apache 2.0 license](https://www.apache.org/licenses/LICENSE-2.0.html).
