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

gsl_histogram *spectrum_read(const char *filename, size_t skip, size_t channels_max, size_t column, size_t compress);
gsl_histogram *spectrum_read_detector(const char *filename, const detector *det);
void spectrum_set_calibration(gsl_histogram *h, const detector *det, int Z);
double spectrum_roi(const gsl_histogram *h, size_t low, size_t high); /* low and high are both inclusive channel numbers*/
size_t spectrum_channels_in_range(const gsl_histogram *h, size_t low, size_t high);
int spectrum_compare(const gsl_histogram *h1, const gsl_histogram *h2, size_t low, size_t high, double *out);
#endif // JABS_SPECTRUM_H
