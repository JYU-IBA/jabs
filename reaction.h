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
#include "ion.h"

typedef struct {
    const jibal_isotope *isotope; /* target isotope */
    int i_isotope; /* location of target isotope in concentration table */
    ion p; /* reaction product, e.g. in case of He-RBS this will be He. */
    double K;
    double max_depth;
    /* TODO: type? */
    /* TODO: cross section to use? */
    /* TODO: cross sections from files */
    double E; /* Previous energy bin, TODO: make a histogram instead! */
    double S; /* Previous straggling bin */
    int stop;
} reaction;

void reactions_print(FILE *f, reaction *reactions, int n_reactions);
#endif //JABS_REACTION_H
