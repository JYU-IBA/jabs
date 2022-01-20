/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021-2022 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#ifndef JABS_APERTURE_H
#define JABS_APERTURE_H

#include <jibal.h>


typedef enum {
    APERTURE_NONE = 0,
    APERTURE_CIRCLE = 1,
    APERTURE_RECTANGLE = 2
} aperture_type;

static const jibal_option aperture_option[] = {
        {JIBAL_OPTION_STR_NONE, APERTURE_NONE},
        {"circle", APERTURE_CIRCLE},
        {"rectangle", APERTURE_RECTANGLE},
        {NULL, 0}
};

typedef struct aperture {
    aperture_type type;
    double diameter;
    double width;
    double height;
} aperture;

const char *aperture_name(const aperture *a);
aperture *aperture_default();
void aperture_free(aperture *a);
double aperture_width_shape_product(const aperture *a, char direction);
aperture *aperture_from_argv(const jibal *jibal, int * const argc, char * const ** const argv);

#endif // JABS_APERTURE_H
