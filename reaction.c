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


void reactions_print(FILE *f, const reaction *reactions) {
    int i = 1;
    for (const reaction *r = reactions; r->type != REACTION_NONE; r++) {
        fprintf(f, "Reaction %3i: %s with %5s.\n", i++, reaction_name(r), r->target->name);
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
        case REACTION_ARB:
            return "ARB";
        default:
            return "???";
    }
}

size_t reaction_count(const reaction *reactions) {
    size_t n=0;
    if(!reactions)
        return 0;
    for (const reaction *r = reactions; r->type != REACTION_NONE; r++) {
        n++;
    }
    return n;
}



reaction reaction_make(const jibal_isotope *incident, const jibal_isotope *target, reaction_type type, jibal_cross_section_type cs, double theta) {
    reaction r;
    r.type = type;
    r.cs = cs;
    r.incident = incident;
    r.target = target;
    if(!r.target) {
        r.type = REACTION_NONE;
        return r;
    }
    if(type == REACTION_RBS) {
        double theta_max=asin(r.target->mass/r.incident->mass);
        if(r.incident->mass >= r.target->mass && theta > theta_max) {
            fprintf(stderr, "RBS with %s is not possible (theta max %g deg, sim theta %g deg)\n", r.target->name, theta_max/C_DEG, theta/C_DEG);
            r.type = REACTION_NONE;
            return r;
        }
    } else if (type == REACTION_ERD) {
        if(theta > C_PI/2.0) {
            fprintf(stderr, "ERD with %s is not possible (theta %g deg > 90.0 deg)", r.target->name, theta);
            r.type = REACTION_NONE;
            return r;
        }
    } else {
        r.type = REACTION_NONE;
        return r;
    }
    return r;
}
