/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#ifndef JABS_SPECTRUM_H
#define JABS_SPECTRUM_H

#include "simulation.h"

gsl_histogram *read_experimental_spectrum(const char *filename, size_t n);
void set_experimental_spectrum_calibration(gsl_histogram *h, const simulation *sim);

#endif // JABS_SPECTRUM_H