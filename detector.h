/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021-2022 Jaakko Julin

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
#include <jibal_masses.h>
#include "sample.h"
#include "aperture.h"
#include "rotate.h"

typedef enum {
    DETECTOR_NONE = 0,
    DETECTOR_ENERGY = 1,
    DETECTOR_TOF = 2,
    DETECTOR_ELECTROSTATIC = 3
} detector_type;

static const jibal_option detector_option[] = {
        {JIBAL_OPTION_STR_NONE, DETECTOR_NONE},
        {"energy",              DETECTOR_ENERGY},
        {"tof",                 DETECTOR_TOF},
        {"electrostatic",       DETECTOR_ELECTROSTATIC},
        {NULL,                  0}
};

typedef struct detector {
    detector_type type;
    double slope;
    double offset;
    double length; /* For ToF */
    double resolution; /* Stored as variance, i.e. energy squared (in SI-units J^2) */
    double theta; /* Polar angle [0, pi] */
    double phi; /* Azimuthal angle [0, 2pi] */
    double solid;
    aperture *aperture;
    double distance;
    size_t column;
    size_t channels;
    size_t compress;
    sample_model *foil_sm;
    sample *foil;
} detector;

inline double detector_calibrated(const detector *det, size_t ch) {return det->offset + det->slope * (unsigned int)(ch*det->compress);}
const char *detector_type_name(const detector *det);
int detector_sanity_check(const detector *det);
detector *detector_from_file(const jibal *jibal, const char *filename); /* one detector from JIBAL configuration style file */
detector **detectors_from_file(const jibal *jibal, const char *filename, size_t *n_detectors_out); /* multiple detectors in a table. TODO: work in progress */
detector *detector_default(detector *det); /* if det is NULL, this returns pointer to a newly allocated det */
void detector_free(detector *det);
int detector_print(const char *filename, const detector *det);
int detector_aperture_set_from_argv(const jibal *jibal, detector *det, int *argc, char * const **argv);
int detector_foil_set_from_argv(const jibal *jibal, detector *det, int *argc, char * const **argv);
int detector_update_foil(detector *det);
int detector_set_var(const jibal *jibal, detector *det, const char *var_str, const char *val_str);
jibal_config_var *detector_make_vars(detector *det);
double detector_angle(const detector *det, char direction);
double detector_theta_deriv(const detector *det, char direction);
double detector_solid_angle_calc(const detector *det);
double detector_resolution(const detector *det, const jibal_isotope *isotope, double E);
#endif //JABS_DETECTOR_H
