/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#ifndef JABS_DETECTOR_H
#define JABS_DETECTOR_H
#include <ctype.h>
#include <jibal_units.h>
#include "sample.h"

typedef struct detector {
    double slope;
    double offset;
    double resolution; /* Stored as variance, i.e. energy squared (in SI-units J^2) */
    double theta; /* Polar angle [0, pi] */
    double phi; /* Azimuthal angle [0, 2pi] */
    double solid;
    size_t number;
    size_t channels;
    size_t compress;
    sample *foil;
    char *foil_description;
} detector;

inline double detector_calibrated(const detector *det, size_t ch) {return det->offset + det->slope * (unsigned int)(ch*det->compress);}
int detector_sanity_check(const detector *det);
detector *detector_from_file(const jibal *jibal, const char *filename);
detector *detector_default();
void detector_free(detector *det);
int detector_print(const char *filename, const detector *det);
int detector_update_foil(const jibal *jibal, detector *det);
int detector_set_var(const jibal *jibal, detector *det, const char *var_str, const char *val_str);
jibal_config_var *detector_make_vars(detector *det);
#endif //JABS_DETECTOR_H
