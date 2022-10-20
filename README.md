# JaBS - Jaakko's Backscattering Simulator

[![DOI](https://zenodo.org/badge/414092526.svg)](https://zenodo.org/badge/latestdoi/414092526)

Simulates and fits RBS, EBS, ERD and NRA spectra rapidly.
    
    Copyright (C) 2021 - 2022 Jaakko Julin
    
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
Some ready-to-use [binary packages](http://users.jyu.fi/~jaakjuli/jabs/binaries/) for Windows and macOS may be available. The Windows distribution includes two executables: jabs.exe (command line) and qjabs.exe (graphical interface).

See instructions below on how to compile JaBS from sources. Supported platforms are Windows, macOS and Linux.

## Getting started

One [example](example/) showing basic use is provided, along with more complex [test-cases](example/tests/). Try running `help` and `help commands` in the interactive mode.

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

The interactive mode should be self-explanatory and an internal help is provided. Please see [the example script](example/example.jbs) to get started. Note that scripts are executed in the directory where the file is located, when opened in the GUI, but in the current working directory by the command-line program. This affects loading and saving filenames.

The scripting language is not a programming language, there is no flow control, new variables can not be introduced etc. This may change in the future.

## Command line usage

Several parameters can be set from the command line. See `jabs -h` and try something like this:

~~~~
$ ./jabs -E 2MeV --alpha=10deg --theta=170deg --out=spectrum.csv sample Au 500tfu SiO2 1000tfu Si 10000tfu
~~~~

Detector and sample can be read from files. The file formats are simple and human readable. Please see [the example](example).

## Features
### Implemented
 - Basic RBS spectrum simulation with Rutherford and Andersen cross-sections.
 - Atomic data, electronic stopping and straggling given by [JIBAL](https://github.com/JYU-IBA/jibal).
 - Automatic recoil spectra when working in forward angles (ERDA)
 - Arbitrary cross sections and reactions from R33 files. Both EBS (Q-value is zero) as well as p-p NRA are implemented. 
 - Point-by-point and layered sample models
 - Roughness (gamma distribution)
 - Arbitrary geometry, detector and sample tilt can be expressed in arbitrary spherical coordinates
 - Non-linear detector calibration (polynomial of arbitrary degree) and constant detector energy resolution or timing resolution (for ToF detectors)
 - Different calibrations are possible for different elements (Z-specific calibration).
 - Multilayer foil in front of detector
 - Reading of experimental data, integer factor compress
 - Fitting to experimental spectrum (multidimensional non-linear least squares fitting)
 - Linear ad-hoc substrate channeling correction
 - Faster (less accurate) modes
 - Slower (more accurate) modes
 - User defined "molecules" i.e. elements with fixed concentration ratios (e.g. you can fit C in SiO2 without changing Si/O ratio)
 - Weighting of cross-sections by straggling 
 - Kinematic (geometric) broadening due to finite detector size and beam spot
 - Dual scattering model (although it is not guaranteed to be accurate)
 - Simultaneous multi-detector simulation and fitting.
 - Stopping corrections can be supplied by user (Bragg correction) for a specific layer as well as similar corrections for straggling and yield
 - Different detector calibration can be given for different particles (different proton number Z)
 - Two-phase fitting, where faster physics model is used in the beginning and more accurate (user configurable) model is turned on after the first phase is starting to converge. User can skip the faster fitting phase.
 - Testing of areal sum (counts) and residuals. Some test cases are run by the developer to check sanity and accuracy of simulations for every new release.
 - Higher accuracy mode using adaptive integration for more accurate handling of sharp peaks in cross sections (resonances) and accurate weighting of cross sections by (Gaussian) straggling.

### Not (yet) implemented, but planned
 - Support for more input and output data formats (CSV, IDF, ...)
 - Multiple scattering (small angle)
 - Simulation of pile-up and dead time
 - Simulation of time-of-flight spectra
 - Publication quality plotting

### Distant future, if ever
 - Parallel processing of independent simulations
 - Advanced fitting algorithms
 - Fitting of spectra from different measurements (different beam, fluence etc for each simulation)
 - Turing completeness of the scripting language
### Known issues
 - Saving a detector to a file is not supported. Calibrations can be saved as a script.
 - Transmission geometry is not supported
 - Ad-hoc channeling correction is the same for all detectors
 - Dual scattering assumes first scattering is RBS (not ERD). Cross sections are not calculated accurately (must use integrated cross sections instead of approximating using differential cross sections since solid angles involved are large).
 - Dual scattering is benchmarked against SimNRA and is known to produce somewhat different results.

## Fitting

There is a fitting feature, activated with the `--fit` or `-F` option or `fit` script command. The multidimensional nonlinear least-squares fitting is based on [GSL multifit](https://www.gnu.org/software/gsl/doc/html/nls.html). The accuracy and sanity of fits must be evaluated by the user.

The command line interface for fitting is currently undocumented and may be removed in the future.

Interactive/script example:

    jabs> add fit range [500:900] [1000:1200]
    jabs> fit *calib*,fluence,thick1
    jabs> show fit

Use `show fit variables` to see a list of variables that can be fitted. This list changes when the sample and detector(s) change. Using `*` or `?` wildcards is possible, the example above will add all calibration parameters of all detectors.

## Open source licenses

This repository contains following open source software, or parts of them:
 * Code derived from GSL (GNU Scientific Library) examples, e.g. in [fit.c](fit.c). Used under the terms of GPL v3.
 * QCustomPlot [qcustomplot.cpp](qjabs/qcustomplot.cpp) and [qcustomplot.h](qjabs/qcustomplot.h), Copyright (C) 2011-2021 Emanuel Eichhammer. Used under the terms of GPL v3.
 * Coordinate rotate routine from MCERD [rotate.c](rotate.c). Copyright Kai Arstila. Used under the terms of GPL v2.
 * Some routines from NetBSD, Copyright (c) 2011 The NetBSD Foundation, Inc. See [win_compat.c](win_compat.c) for the full copyright notice.
 * [Qt6](https://www.qt.io/) GUI is distributed under the provisions of [LGPL version 3](https://doc.qt.io/qt-6/lgpl.html). Qt Toolkit is Copyright (C) 2018 The Qt Company Ltd. and other contributors.
 * Qt Toolkit examples used under terms of the BSD license
 * [Material Design](https://google.github.io/material-design-icons/) [icons](qjabs/icons) licensed under [Apache 2.0 license](https://www.apache.org/licenses/LICENSE-2.0.html).
 * Version tracking from Git, [cmake-git-version-tracking](https://github.com/andrew-hardin/cmake-git-version-tracking) Copyright (c) 2020 Andrew Harding. Used under the terms of [MIT License](https://raw.githubusercontent.com/andrew-hardin/cmake-git-version-tracking/master/LICENSE).
