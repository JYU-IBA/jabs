/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2022 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#ifndef CALIB_CALIBRATION_H
#define CALIB_CALIBRATION_H
#include <inttypes.h>
#include <stdlib.h>
#include <jibal_option.h>
#include "histogram.h"

typedef enum calibration_type {
    CALIBRATION_NONE = 0,
    CALIBRATION_ARB = 1,
    CALIBRATION_LINEAR = 2,
    CALIBRATION_POLY = 3
} calibration_type;

static const jibal_option calibration_option[] = {
        {JIBAL_OPTION_STR_NONE, CALIBRATION_NONE},
        {"arbitrary", CALIBRATION_ARB}, /* Not implemented */
        {"linear", CALIBRATION_LINEAR},
        {"poly", CALIBRATION_POLY},
        {NULL, 0}
};

typedef enum calibration_param_type {
    CALIBRATION_PARAM_RESOLUTION = -1,
    CALIBRATION_PARAM_OFFSET = 0,
    CALIBRATION_PARAM_SLOPE = 1,
    CALIBRATION_PARAM_QUAD = 2 /* Higher order terms etc are possible */
} calibration_param_type;

typedef struct calibration {
    calibration_type type;
    double (*f)(const void *, size_t);
    void *params;
    double resolution; /* Stored as FWHM in relevant SI units. Note that can be e.g. energy or time depending on detector type. */
    double resolution_variance; /* Calculated based on "resolution" before needed by detector_update(). */
} calibration;

typedef struct calibration_params_linear {
    double offset;
    double slope;
} calibration_params_linear;

typedef struct calibration_params_poly {
    size_t n;
    double *a; /* Array, size n+1, polynomial is f(x) = a[0] + a[1] * x + a[2] * x^2 + ... + a[n] * x^n */
} calibration_params_poly;

calibration *calibration_init(void);
void calibration_free(calibration *c);
calibration *calibration_clone(const calibration *c_orig);
calibration *calibration_init_linear(void);
calibration *calibration_init_poly(size_t n); /* n is the degree of the polynomial */
double calibration_linear(const void *params, size_t x);
double calibration_poly(const void *params, size_t x);
double calibration_none(const void *params, size_t x);
inline double calibration_eval(const calibration *c, size_t ch) {return c->f(c->params, ch);}
size_t calibration_inverse(const calibration *cal, double E, size_t ch_max); /* Returns number between [0,ch_max), returns 0 also if E is outside calibration. calibration curve has to be monotonous */
int calibration_set_param(calibration *c, int i, double value);
size_t calibration_get_number_of_params(const calibration *c);
double calibration_get_param(const calibration *c, int i); /* get i'th param in range [0..n-1], get n by  calibration_get_number_of_params()*/
double *calibration_get_param_ref(calibration *c, int i);
int calibration_copy_params(calibration *dst, const calibration *src); /* Copies parameters from src to dst. Calibration types may be different (e.g. linear or poly), we'll do our best. */
const char *calibration_name(const calibration *c);
char *calibration_to_string(const calibration *c);
char *calibration_param_name(calibration_type type, calibration_param_type i); /* Name of i'th param (e.g. "slope"). Returns a string that can be free'd */
int calibration_is_monotonically_increasing(const calibration *cal, size_t n_channels);
void calibration_apply_to_histogram(const calibration *cal, jabs_histogram *h);
#endif //CALIB_CALIBRATION_H
