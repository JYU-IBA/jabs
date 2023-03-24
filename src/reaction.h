/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_REACTION_H
#define JABS_REACTION_H

#include <jibal_masses.h>
#include <jibal.h>
#include <jibal_cross_section.h>
#include <jibal_kin.h>
#include <jibal_r33.h>
#include "ion.h"
#include "plugin.h"
#include "sim_calc_params.h"

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

typedef enum jabs_reaction_cs {
    JABS_CS_NONE = 0,
    JABS_CS_RUTHERFORD = 1,
    JABS_CS_ANDERSEN = 2,
    JABS_CS_UNIVERSAL = 3
} jabs_reaction_cs;

static const jibal_option jabs_cs_types[] = {
        {JIBAL_OPTION_STR_NONE, JABS_CS_NONE},
        {"Rutherford", JABS_CS_RUTHERFORD},
        {"Andersen", JABS_CS_ANDERSEN},
        {"Universal", JABS_CS_UNIVERSAL},
        {NULL, 0}
};

typedef struct reaction {
    reaction_type type;
    char *name;
    /* Reactions are like this: target(incident,product)residual */
    const jibal_isotope *incident;
    const jibal_isotope *target;
    const jibal_isotope *product;
    const jibal_isotope *residual; /* Typically the "heavy" reaction product (of limited interest) */
    jabs_reaction_cs cs; /* Cross section model to use (e.g. screening corrections) */
    char *filename; /* for REACTION_FILE and REACTION_PLUGIN */
#ifdef JABS_PLUGINS
    jabs_plugin *plugin; /* for REACTION_PLUGIN */
    jabs_plugin_reaction *plugin_r; /* for REACTION_PLUGIN */
#endif
    struct reaction_point *cs_table; /* for REACTION_FILE */
    size_t n_cs_table;
    double theta; /* For REACTION_FILE */
    double Q;
    double E_min;
    double E_max;
    nuclear_stopping *nucl_stop; /* Nuclear stopping for reaction product, may be shared. Can be calculated for a reaction, since this is quite universal. */
    jabs_ion_gsto *ion_gsto;
} reaction;


void reactions_print(reaction * const *reactions, size_t n_reactions);
reaction *reaction_make(const jibal_isotope *incident, const jibal_isotope *target, reaction_type type, jabs_reaction_cs cs);
reaction *reaction_make_from_argv(const jibal *jibal, const jibal_isotope *incident, int *argc, char * const **argv);
const char *reaction_name(const reaction *r);
const char *reaction_type_to_string(reaction_type type);
reaction_type reaction_type_from_string(const char *s);
void reaction_free(reaction *r);
int reaction_is_possible(const reaction *r, const sim_calc_params *p, double theta); /* If this returns FALSE, reaction r is not possible with lab angle theta. If it returns TRUE, it might be. */
int reaction_theta_within_tolerance(const reaction *r, const sim_calc_params *p, double theta);
reaction *r33_file_to_reaction(const jibal_isotope *isotopes, const r33_file *rfile);
int reaction_compare(const void *a, const void *b);
double reaction_product_energy(const reaction *r, double theta, double E); /* Hint: call with E == 1.0 to get kinematic factor */
const char *jabs_reaction_cs_to_string(jabs_reaction_cs cs);
int reaction_generate_name(reaction *r); /* Used internally, don't call. */
jabs_reaction_cs jabs_reaction_cs_from_jibal_cs(jibal_cross_section_type jcs);
#endif //JABS_REACTION_H
