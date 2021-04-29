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
    REACTION_FILE = 3,
    REACTION_ARB = 4 /* TODO: types of reactions */
} reaction_type;

typedef struct reaction {
    reaction_type type;
    const jibal_isotope *incident;
    const jibal_isotope *target; /* target isotope */
    jibal_cross_section_type cs;
    char *filename; /* for REACTION_FILE */
} reaction;


void reactions_print(FILE *f, const reaction *reactions);
reaction reaction_make(const jibal_isotope *incident, const jibal_isotope *target, reaction_type type, jibal_cross_section_type cs, double theta); /* theta is used to check sanity of RBS and ERD reactions */
const char *reaction_name(const reaction *r);
size_t reaction_count(const reaction *reactions);
#endif //JABS_REACTION_H
