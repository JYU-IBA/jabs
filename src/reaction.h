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
#include <jibal.h>
#include <jibal_cross_section.h>
#include <jibal_kin.h>
#include <jibal_r33.h>
#include "ion.h"
#include "plugin.h"

typedef enum {
    REACTION_NONE = 0,
    REACTION_RBS = 1,
    REACTION_RBS_ALT = 2,
    REACTION_ERD = 3,
    REACTION_FILE = 4,
    REACTION_PLUGIN = 5
} reaction_type;

struct reaction_point {
    double E;
    double sigma;
};

typedef struct reaction {
    reaction_type type;
    /* Reactions are like this: target(incident,product)product_nucleus */
    const jibal_isotope *incident;
    const jibal_isotope *target;
    const jibal_isotope *product;
    const jibal_isotope *product_nucleus; /* Not used */
    jibal_cross_section_type cs; /* Cross section model to use (e.g. screening corrections) */
    char *filename; /* for REACTION_FILE and REACTION_PLUGIN */
    struct reaction_point *cs_table; /* for REACTION_FILE */
    jabs_plugin *plugin; /* for REACTION_PLUGIN */
    jabs_plugin_reaction *plugin_r; /* for REACTION_PLUGIN */
    size_t n_cs_table;
    double theta; /* For REACTION_FILE */
    double Q;
    double E_min;
    double E_max;
} reaction;


void reactions_print(FILE *f, reaction * const *reactions, size_t n_reactions);
reaction *reaction_make(const jibal_isotope *incident, const jibal_isotope *target, reaction_type type, jibal_cross_section_type cs);
reaction *reaction_make_from_argv(const jibal *jibal, const jibal_isotope *incident, int *argc, char * const **argv);
const char *reaction_name(const reaction *r);
reaction_type reaction_type_from_string(const char *s);
void reaction_free(reaction *r);
int reaction_is_same(const reaction *r1, const reaction *r2); /* TRUE (=1) if r1 and r2 describe the same reaction. Note that "type" can be different. */
reaction *r33_file_to_reaction(const jibal_isotope *isotopes, const r33_file *rfile);
reaction *plugin_reaction(jabs_plugin *plugin);
int reaction_compare(const void *a, const void *b);
double reaction_product_energy(const reaction *r, double theta, double E); /* Hint: call with E == 1.0 to get kinematic factor */
#endif //JABS_REACTION_H
