/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include "reaction.h"

void reactions_print(FILE *f, const reaction *reactions) {
    int i = 1;
    for (const reaction *r = reactions; r->type != REACTION_NONE; r++) {
        fprintf(f, "Reaction %3i: %s with %5s (isotope id %3i): K = %7.5lf, max depth = %9.3lf tfu\n",
                i++, reaction_name(r), r->isotope->name, r->i_isotope, r->K, r->max_depth / C_TFU);
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
