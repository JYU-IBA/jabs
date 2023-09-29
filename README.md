# JaBS - Jaakko's Backscattering Simulator

[![DOI](https://zenodo.org/badge/414092526.svg)](https://zenodo.org/badge/latestdoi/414092526)

Simulates and fits RBS, EBS, ERD and NRA spectra rapidly.
    
    Copyright (C) 2021 - 2023 Jaakko Julin
    
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
Some ready-to-use binary packages for each [release](https://github.com/JYU-IBA/jabs/releases) for Windows and macOS may be available. The Windows distribution includes two executables: `jabs.exe` (command line) and `qjabs.exe` (graphical interface). The macOS binary is a regular application bundle. All necessary data files (stopping etc.) are included in these distributions.

See instructions below on how to compile JaBS from sources. Supported platforms are Windows, macOS and Linux.

## Getting started

One [example](example/) showing basic use is provided, along with more complex [test-cases](example/tests/). Try running `help` and `help commands` in the interactive mode. XML files in IDF format created by other programs can be imported.

## Installation from sources

See [instructions](INSTALL.md).

## Graphical user interface

Simple Qt-based [GUI](qjabs/) is provided. Scripts can be edited and run, simulated spectra and experimental spectra are plotted (if available).
      
## Interactive or scripted usage

Launch JaBS in interactive mode simply by running `jabs` or `jabs --interactive`. Scripts (i.e. files with commands identical to interactive mode input) can be given as command line parameters e.g. `jabs script.jbs`, piped in `jabs < script.jbs` or run using the interactive mode.

The interactive mode should be self-explanatory and an internal help is provided. Please see [the example script](example/example.jbs) to get started. Note that scripts are executed in the directory where the file is located, when opened in the GUI, but in the current working directory by the command-line program. This affects loading and saving filenames.

The scripting language is not a programming language, there is no flow control, new variables can not be introduced etc. However multiple simulations/fits are possible in a single file and other scripts may be loaded, which in turn may load other scripts.

## Features
### Implemented
 - Basic RBS spectrum simulation with Rutherford and Andersen cross-sections.
 - Atomic data, electronic stopping and straggling given by [JIBAL](https://github.com/JYU-IBA/jibal).
 - Automatic simulation of recoil spectra when working in forward angles (ERDA) and simulation of both solutions (+/-) for RBS when possible
 - Arbitrary cross sections and reactions from R33 files. Both EBS (Q-value is zero) and p-p NRA are implemented. 
 - Point-by-point and layered sample models
 - Roughness using a gamma distribution and arbitrary roughness using files (weight and thickness tables)
 - Arbitrary geometry, detector and sample tilt can be expressed in arbitrary spherical coordinates
 - Transmission geometry
 - Non-linear detector calibration (polynomial of arbitrary degree) and constant detector energy resolution or timing resolution (for ToF detectors)
 - Different calibrations are possible for different elements (Z-specific calibration).
 - Multilayer foil in front of detector
 - Reading of experimental data, integer factor compress
 - Fitting to experimental spectrum (multidimensional non-linear least squares fitting)
 - Simultaneous multi-detector simulation and fitting.
 - Faster (less accurate) modes
 - Slower (more accurate) modes 
 - Two-phase fitting, where faster physics model is used in the beginning and more accurate (user configurable) model is turned on after the first phase is starting to converge. User can skip the faster fitting phase.
 - User defined "molecules" i.e. elements with fixed concentration ratios (e.g. you can fit C in SiO2 without changing Si/O ratio)
 - Weighting of cross-sections by straggling 
 - Kinematic (geometric) broadening due to finite detector size and beam spot
 - Stopping corrections can be supplied by user (Bragg correction) for a specific layer as well as similar corrections for straggling and yield
 - Testing of areal sum (counts) and residuals. Some test cases are run by the developer to check sanity and accuracy of simulations for every new release.
 - Higher accuracy mode using adaptive integration for more accurate handling of sharp peaks in cross sections (resonances) and accurate weighting of cross sections by (Gaussian) straggling.
 - Conversion tool from IDF to JaBS script (partial support)
 - Simulation of large angle plural scattering (dual scattering model), with the assumption that first scattering is RBS (not ERD).
 - Multiprocessor support when simulating multiple detectors and when fitting multiple parameters
 - 
### Not implemented, but planned or being worked on
 - Support for more input and output data formats (CSV, ...)
 - Multiple scattering (small angle)
 - Simulation of pile-up and dead time
 - Simulation of time-of-flight spectra
 - Publication quality plotting

### Distant future, if ever
 - Advanced fitting algorithms
 - Fitting of spectra from different measurements (different beam, fluence etc for each simulation)
 - Turing completeness of the scripting language

### Known issues
 - Saving a detector to a file is not supported. Calibrations can be saved as a script.
 - Ad-hoc channeling correction is the same for all detectors
 - Import of IDF is partial at best and exporting is not implemented

## Fitting

The multidimensional nonlinear least-squares fitting is based on Levenberg-Marquardt algorithm [GSL multifit](https://www.gnu.org/software/gsl/doc/html/nls.html).

Example:

    jabs> add fit range [500:900] [1000:1200]
    jabs> fit *calib*,fluence,thick1
    jabs> show fit

Use `show fit variables` to see a list of variables that can be fitted. This list changes when the sample and detector(s) change. Using `*` or `?` wildcards is possible, the example above will add all calibration parameters of all detectors.

## Open source licenses

This repository contains following open source software, or parts of them:
 * Code derived from GSL (GNU Scientific Library) examples, e.g. in [fit.c](src/fit.c). Used under the terms of GPL v3.
 * QCustomPlot [qcustomplot.cpp](qjabs/qcustomplot.cpp) and [qcustomplot.h](qjabs/qcustomplot.h), Copyright (C) 2011-2021 Emanuel Eichhammer. Used under the terms of GPL v3.
 * Coordinate rotate routine from MCERD [rotate.c](src/rotate.c). Copyright Kai Arstila. Used under the terms of GPL v2.
 * Some routines from NetBSD, Copyright (c) 2011 The NetBSD Foundation, Inc. See [win_compat.c](src/win_compat.c) for the full copyright notice.
 * Qt Toolkit examples used under terms of the BSD license.
 * [Material Design](https://google.github.io/material-design-icons/) [icons](qjabs/icons) licensed under [Apache 2.0 license](https://www.apache.org/licenses/LICENSE-2.0.html).
 * Version tracking from Git, [cmake-git-version-tracking](https://github.com/andrew-hardin/cmake-git-version-tracking) Copyright (c) 2020 Andrew Harding. Used under the terms of [MIT License](https://raw.githubusercontent.com/andrew-hardin/cmake-git-version-tracking/master/LICENSE).

Binary distributions contain following open source software:
* [Qt6](https://www.qt.io/) GUI is distributed under the provisions of [LGPL version 3](https://doc.qt.io/qt-6/lgpl.html). Qt Toolkit is Copyright (C) 2018 The Qt Company Ltd. and other contributors.
 * [libxml2](https://gitlab.gnome.org/GNOME/libxml2) library, used under the terms of [MIT License](https://www.opensource.org/licenses/mit-license.html). Copyright (C) 1998-2012 Daniel Veillard. 
