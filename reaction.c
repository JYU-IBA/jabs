/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include <assert.h>
#include <string.h>
#include "defaults.h"
#include "reaction.h"
#include "message.h"


void reactions_print(FILE *f, reaction * const * reactions, size_t n_reactions) {
    if(!reactions)
        return;
    jabs_message(MSG_INFO, f, "Reactions (%zu):\n", n_reactions);
    for(size_t i = 0; i < n_reactions; i++) {
        const reaction *r = reactions[i];
        if(!r) {
#ifdef DEBUG
            fprintf(stderr, "Reaction i=%zu (%zu) is NULL\n", i, i + 1);
#endif
            continue;
        }
        jabs_message(MSG_INFO, f, "%3zu: %4s with %5s (reaction product %s).", i + 1, reaction_name(r), r->target->name, r->product->name);
        if(r->type == REACTION_FILE) {
            jabs_message(MSG_INFO, f, " Incident = %s, Theta = %g deg, E = [%g keV, %g keV]. Q = %g MeV. Data from file \"%s\".\n", r->incident->name, r->theta/C_DEG, r->cs_table[0].E/C_KEV, r->cs_table[r->n_cs_table-1].E/C_KEV, r->Q/C_MEV, r->filename);
        } else {
            jabs_message(MSG_INFO, f, " %s cross sections (built-in).", jibal_cs_types[r->cs]);
            if(r->E_max < E_MAX || r->E_min > 0.0) {
                jabs_message(MSG_INFO, f, " E = [%g keV, %g keV]", r->E_min/C_KEV, r->E_max/C_KEV);
            }
            jabs_message(MSG_INFO, f, "\n");
        }
    }
}

const char *reaction_name(const reaction *r) {
    if(!r)
        return "NULL";
    switch (r->type) {
        case REACTION_NONE:
            return "NONE";
        case REACTION_RBS:
            return "RBS";
        case REACTION_ERD:
            return "ERD";
        case REACTION_FILE:
            return "FILE";
        case REACTION_ARB:
            return "ARB";
        default:
            return "???";
    }
}

reaction_type reaction_type_from_string(const char *s) {
    if(!s)
        return REACTION_NONE;
    if(strcmp(s, "RBS") == 0)
        return REACTION_RBS;
    if(strcmp(s, "ERD") == 0)
        return REACTION_ERD;
    if(strcmp(s, "FILE") == 0)
        return REACTION_FILE;
#ifdef DEBUG
    fprintf(stderr, "Reaction type string %s is not valid.\n", s);
#endif
    return REACTION_NONE;
}

reaction *reaction_make(const jibal_isotope *incident, const jibal_isotope *target, reaction_type type, jibal_cross_section_type cs) {
    if(!target || !incident) {
        return NULL;
    }
#if 0
    if(!force) {
        if(type == REACTION_RBS) {
            double theta_max = asin(target->mass / incident->mass);
            if(incident->mass >= target->mass && theta > theta_max) {
                fprintf(stderr, "RBS with %s is not possible (theta max %g deg, sim theta %g deg)\n", target->name,
                        theta_max / C_DEG, theta / C_DEG);
                return NULL;
            }
        } else if(type == REACTION_ERD) {
            if(theta > C_PI / 2.0) {
                fprintf(stderr, "ERD with %s is not possible (theta %g deg > 90.0 deg)", target->name, theta);
                return NULL;
            }
        }
    }
#endif
    reaction *r = malloc(sizeof(reaction));
    r->type = type;
    r->cs = cs;
    r->incident = incident;
    r->target = target;
    r->filename = NULL;
    r->cs_table = NULL;
    r->n_cs_table = 0;
    r->E_min = 0.0;
    r->E_max = E_MAX;
    r->Q = 0.0;
    switch(type) {
        case REACTION_ERD:
            r->product = target;
            r->product_nucleus = incident;
            break;
        case REACTION_RBS:
            r->product = incident;
            r->product_nucleus = target;
            break;
        default:
            fprintf(stderr, "Warning, reaction product is null.\n");
            r->product = NULL;
            break;
    }
    return r;
}

reaction *reaction_make_from_argv(const jibal *jibal, const jibal_isotope *incident, int argc, char * const *argv) {
    if(argc < 2) {
        jabs_message(MSG_ERROR, stderr, "Not enough arguments\n");
        return NULL;
    }
    reaction_type type = reaction_type_from_string(argv[0]);
    const jibal_isotope *target = jibal_isotope_find(jibal->isotopes, argv[1], 0, 0);
    if(type == REACTION_NONE) {
        jabs_message(MSG_ERROR, stderr, "This is not a valid reaction type: \"%s\".\n", argv[0]);
        return NULL;
    }
    if(!target) {
        jabs_message(MSG_ERROR, stderr, "This is not a valid isotope: \"%s\".\n", argv[1]);
        return NULL;
    }
    reaction *r = reaction_make(incident, target, type, JIBAL_CS_NONE); /* Warning: JIBAL_CS_NONE used here, something sane must be supplied after this somewhere! */
    argc -= 2;
    argv += 2;
    while (argc >= 2) {
        if(strcmp(argv[0], "max") == 0) {
            r->E_max = jibal_get_val(jibal->units, 'E', argv[1]);
        } else if(strcmp(argv[0], "min") == 0) {
            r->E_min = jibal_get_val(jibal->units, 'E', argv[1]);
        } else if(strcmp(argv[0], "cs") == 0) {
            r->cs =  jibal_option_get_value(jibal_cs_types, argv[1]);
        }
        argc -= 2;
        argv += 2;
    }
    return r;
}

void reaction_free(reaction *r) {
    if(!r)
        return;
    free(r->cs_table);
}

int reaction_is_same(const reaction *r1, const reaction *r2) {
    if(!r1 || !r2)
        return FALSE;
    if(r1->incident != r2->incident)
        return FALSE;
    if(r1->target != r2->target)
        return FALSE;
    if(r1->product != r2->product)
        return FALSE;
#if 0 /* TODO: support for multiple detectors and multiple reactions from files? */
    if(fabs(r1->theta - r2->theta) > 0.01*C_DEG)
        return FALSE;
#endif
    return TRUE;
}

reaction *r33_file_to_reaction(const jibal_isotope *isotopes, const r33_file *rfile) {
    const jibal_isotope *nuclei[R33_N_NUCLEI];
    for(size_t i = 0; i < R33_N_NUCLEI; i++) {
#ifdef R33_IGNORE_REACTION_STRING
        nuclei[i] = jibal_isotope_find(isotopes, NULL, r33_double_to_int(rfile->zeds[i]), r33_double_to_int(rfile->masses[i]));
#else
        nuclei[i] = jibal_isotope_find(isotopes, rfile->reaction_nuclei[i], 0, 0);
#endif
        if(!nuclei[i]) {
#ifdef R33_IGNORE_REACTION_STRING
            fprintf(stderr, "Could not parse an isotope from Z=%g, mass=%g.", rfile->zeds[i], rfile->masses[i]);
#else
            jabs_message(MSG_ERROR, stderr, "Could not parse an isotope from \"%s\".\n", rfile->reaction_nuclei[i]);
#endif
            return NULL;
        }
    }
    if(rfile->composition) {
        jabs_message(MSG_ERROR, stderr, "This program does not currently support \"Composition\" in R33 files.\n");
        return NULL;
    }


    reaction *r = malloc(sizeof(reaction));
    r->cs = JIBAL_CS_NONE;
    r->type = REACTION_FILE;
#ifdef R33_IGNORE_REACTION_STRING
    r->incident = nuclei[0];
    r->target = nuclei[1];
#else
    r->incident = nuclei[1];
    r->target = nuclei[0];
#endif
    r->product = nuclei[2];
    r->product_nucleus = nuclei[3];
    if(rfile->unit == R33_UNIT_RR && r->incident != r->product) {
        jabs_message(MSG_ERROR, stderr, "R33 file is in units of ratio to Rutherford, but reaction product (%s) is not the same as target (%s). I don't know what to do.\n", r->product->name, r->target->name);
        return NULL;
    }
    r->theta = rfile->theta * C_DEG;
    r->Q = rfile->Qvalues[0] * C_KEV; /* TODO: other values? */
    if(rfile->filename) {
        r->filename = strdup(rfile->filename);
    } else {
        r->filename = NULL;
    }
    /* TODO: copy and convert data (check units etc) */
    r->n_cs_table = rfile->n_data;
    r->cs_table = malloc(sizeof(struct reaction_point) * r->n_cs_table);
    for(size_t i = 0; i < rfile->n_data; i++) {
        struct reaction_point *rp = &r->cs_table[i];
        rp->E = rfile->data[i][0] * rfile->enfactors[0] * C_KEV; /* TODO: other factors? */
        if(rfile->unit == R33_UNIT_RR) {
            rp->sigma = rfile->data[i][2] * rfile->sigfactors[0] * jibal_cross_section_rbs(r->incident, r->target, r->theta, rp->E, JIBAL_CS_RUTHERFORD);
        } else if (rfile->unit == R33_UNIT_MB){
            rp->sigma = rfile->data[i][2] * rfile->sigfactors[0] * C_MB_SR;
        } else {
            rp->sigma = 0.0;
        }
    }
    return r;
}

int reaction_compare(const void *a, const void *b) {
    const reaction *r_a = *((const reaction **)a);
    const reaction *r_b = *((const reaction **)b);
    if(r_a == NULL)
        return 1; /* TODO: or 1? */
    if(r_b == NULL)
        return -1; /* TODO: or -1? */
    const jibal_isotope *isotope_a = r_a->target;
    const jibal_isotope *isotope_b = r_b->target;
    if(r_a->type < r_b->type)
        return -1;
    if(r_a->type > r_b->type)
        return 1;
    if(isotope_a->Z == isotope_b->Z) { /* Same Z, compare by A */
        return isotope_a->A - isotope_b->A;
    } else {
        return isotope_a->Z - isotope_b->Z;
    }
}

double reaction_product_energy(const reaction *r, double theta, double E) { /* Hint: call with E == 1.0 to get kinematic factor */
    const double Q = r->Q;
    if(Q == 0.0) {
        if(r->product == r->incident) {
            return jibal_kin_rbs(r->incident->mass, r->target->mass, theta, '+')*E;
        }  else if(r->product == r->target) {
            return jibal_kin_erd(r->incident->mass, r->target->mass, theta)*E;
        } else {
            jabs_message(MSG_WARNING, stderr, "Warning: Reaction Q value is zero for reaction \"%s\" but based on incident %s, target %s and product %s I don't know what to do.\n", r->incident->name, r->target->name, r->product->name);
            return 0.0;
        }
    }
    const double m1 = r->incident->mass;
    const double m2 = r->target->mass;
    const double m3 = r->product->mass;
    const double m4 = r->product_nucleus->mass;
    double E_total = E + Q;
    double a13 = ((m1 * m3)/((m1 + m2)*(m3 + m4))) * (E / E_total);
    double a24 = ((m2 * m4)/((m1 + m2)*(m3 + m4))) * (1 + (m1/m2) * Q / E_total);
    if(a13 > a24) {
        double theta_max = sqrt(asin(a24/a13));
        if(theta > theta_max) {
            return 0.0;
        }
    }
    double E_out = E_total * a13 * pow2(cos(theta) + sqrt(a24/a13 - pow2(sin(theta)))); /* Solution with '-' is ignored. */
    return E_out;
}
