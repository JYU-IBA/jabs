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

void reactions_print(FILE *f, reaction *reactions, int n_reactions) {
    int i;
    for (i = 0; i < n_reactions; i++) {
        reaction *r = &reactions[i];
        fprintf(stderr, "Reaction %3i/%i: RBS with %5s (isotope id %3i): K = %7.5lf, max depth = %9.3lf tfu\n", i+1, n_reactions,
                r->isotope->name, r->i_isotope,
                r->K, r->max_depth / C_TFU);
    }

}
