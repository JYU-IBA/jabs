/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include "ion.h"

void ion_set_isotope(ion *ion, const jibal_isotope *isotope) {
    if(!isotope)
        return;
    ion->isotope = isotope;
    ion->mass = isotope->mass;
    ion->Z = isotope->Z;
}

void ion_set_angle(ion *ion, double angle) {
    if(!ion)
        return;
    ion->angle = angle;
    ion->cosine = cos(angle);
    ion->inverse_cosine = 1.0/ion->cosine;
}

