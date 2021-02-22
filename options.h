/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_OPTIONS_H
#define JABS_OPTIONS_H

#include <jibal.h>

typedef struct {
    jibal *jibal;
    int verbose;
    int fit;
    int fit_low;
    int fit_high;
    char *out_filename;
    char *exp_filename;
    char *bricks_filename;
    char *fit_vars;
} global_options;

#include "simulation.h"

void read_options(global_options *global, simulation *sim, int *argc, char ***argv);
const char *jabs_version();
void usage();

#endif // JABS_OPTIONS_H
