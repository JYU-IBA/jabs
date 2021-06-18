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
#include "reaction.h"


void reactions_print(FILE *f, const reaction *reactions, size_t n_reactions) {
    if(!reactions)
        return;
    for(size_t i = 0; i < n_reactions; i++) {
        const reaction *r = &reactions[i];
        fprintf(f, "Reaction %3zu: %s with %5s (reaction product %s).\n", i + 1, reaction_name(r), r->target->name, r->product->name);
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
    switch(type) {
        case REACTION_ERD:
            r->product = target;
            break;
        case REACTION_RBS:
            r->product = incident;
            break;
        default:
            fprintf(stderr, "Warning, reaction product is null.\n");
            r->product = NULL;
            break;
    }
    return r;
}

void reaction_free(reaction *r) {
    if(!r)
        return;
    free(r->cs_table);
}

int reaction_is_same(const reaction *r1, const reaction *r2) {
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
            fprintf(stderr, "Could not parse an isotope from \"%s\".\n", rfile->reaction_nuclei[i]);
#endif
            return NULL;
        }
    }
    if(rfile->composition) {
        fprintf(stderr, "This program does not currently support \"Composition\" in R33 files.\n");
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
        fprintf(stderr, "R33 file is in units of ratio to Rutherford, but reaction product (%s) is not the same as target (%s). I don't know what to do.\n", r->product->name, r->target->name);
        return NULL;
    }
    r->theta = rfile->theta * C_DEG;
    r->Q = rfile->Qvalues[0]; /* TODO: other values? */
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
