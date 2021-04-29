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
#include "reaction.h"


void reactions_print(FILE *f, reaction * const *reactions) {
    int i = 1;
    if(!reactions)
        return;
    for (reaction * const *r = reactions; *r != NULL; r++) {
        fprintf(f, "Reaction %3i: %s with %5s (reaction product %s).\n", i++, reaction_name(*r), (*r)->target->name, (*r)->product->name);
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

size_t reaction_count(reaction * const *reactions) {
    size_t n=0;
    if(!reactions)
        return 0;
    for (reaction * const *r = reactions; *r != NULL; r++) {
        n++;
    }
    return n;
}



reaction *reaction_make(const jibal_isotope *incident, const jibal_isotope *target, reaction_type type, jibal_cross_section_type cs, double theta, int force) {
    if(!target || !incident) {
        return NULL;
    }
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
