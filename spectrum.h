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

#include "detector.h"
#include "simulation.h"

gsl_histogram *spectrum_read(const char *filename, const detector *det);
void spectrum_set_calibration(gsl_histogram *h, const detector *det, int Z);
double spectrum_roi(gsl_histogram *h, size_t low, size_t high); /* low and high are both inclusive channel numbers*/
size_t spectrum_channels_in_range(gsl_histogram *h, size_t low, size_t high);
#endif // JABS_SPECTRUM_H
