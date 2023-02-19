/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2022 Jaakko Julin

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
#include "calibration.h"

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
    char *name;
    detector_type type;
    struct calibration *calibration;
    struct calibration **calibration_Z; /* Array of calibrations for a given Z. calibration_Z[i] NULL means use default calibration for Z == i. Array is dynamically allocated (number of elements: cal_Z_max + 1) as and when required, so this should be NULL when cal_Z_max == 0. */
    int cal_Z_max;
    double length; /* For ToF */
    double theta; /* Polar angle [0, pi], usually the same as scattering angle (when DS is disabled) */
    double phi; /* Azimuthal angle [0, 2pi] */
    double beta; /* Exit angle, used only if beta calculation is manual (doesn't work with dual scattering or geometric straggling calculation) */
    double solid;
    aperture *aperture;
    double distance;
    size_t column;
    size_t channels;
    size_t compress;
    sample_model *foil_sm;
    sample *foil;
} detector;

calibration *detector_get_calibration(const detector *det, int Z); /* Returns Z specific calibration (if it exists) det->calibration otherwise. */
inline double detector_calibrated(const detector *det, int Z, size_t ch) {return calibration_eval(detector_get_calibration(det, Z), ch*det->compress);}
int detector_set_calibration_Z(const jibal_config *jibal_config, detector *det, calibration *cal, int Z); /* Sets Z specific calibration, (re)allocates space for det->calibration_Z if required. */
const char *detector_type_name(const detector *det);
int detector_sanity_check(const detector *det, size_t n_channels);
detector *detector_default(detector *det); /* if det is NULL, this returns pointer to a newly allocated det */
void detector_free(detector *det);
void detector_calibrations_free(detector *det);
int detector_print(const jibal *jibal, const detector *det);
int detector_aperture_set_from_argv(const jibal *jibal, detector *det, int *argc, char * const **argv);
int detector_foil_set_from_argv(const jibal *jibal, detector *det, int *argc, char * const **argv);
int detector_update_foil(detector *det);
double detector_angle(const detector *det, char direction);
double detector_solid_angle_calc(const detector *det);
double detector_resolution(const detector *det, const jibal_isotope *isotope, double E); /* Isotope is used for Z (Z specific resolution) and for mass (ToF detector) */
char *detector_resolution_to_string(const detector *det, int Z);
void detector_update(detector *det);
const char *detector_param_unit(const detector *det); /* return a suitable unit based on detector type, e.g. "keV" when type == DETECTOR_ENERGY. TODO: use this wisely (we can simulate energy spectra with a ToF detector!) */
double detector_param_unit_factor(const detector *det);
int detector_set_name(detector *det, const char *name);
const char *detector_name(detector *det);
detector *detector_clone(const detector *det_orig);
#endif //JABS_DETECTOR_H
