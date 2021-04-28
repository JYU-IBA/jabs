/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_REACTION_H
#define JABS_REACTION_H

#include <jibal_masses.h>
#include <gsl/gsl_histogram.h>
#include <jibal_cross_section.h>
#include <jibal_kin.h>
#include "ion.h"

typedef enum {
    REACTION_NONE = 0,
    REACTION_RBS = 1,
    REACTION_ERD = 2,
    REACTION_ARB = 3 /* TODO: types of reactions */
} reaction_type;

typedef struct reaction {
    const jibal_isotope *incident;
    const jibal_isotope *target; /* target isotope */
    double K;
    jibal_cross_section_type cs;
    double (*cross_section)(const struct reaction *r);
    /* TODO: cross section to use? */
    /* TODO: cross sections from files */
    reaction_type type;
} reaction;

void reactions_print(FILE *f, const reaction *reactions);
struct reaction reaction_make(const jibal_isotope *incident, const jibal_isotope *target, reaction_type type, double theta);
const char *reaction_name(const reaction *r);
size_t reaction_count(const reaction *reactions);
double reaction_cross_section(const reaction *r);
#endif //JABS_REACTION_H
