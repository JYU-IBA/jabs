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
    int rbs;
    int erd;
    int fit;
    int fit_low;
    int fit_high;
    int print_isotopes;
    int print_iters;
    int ds;
    int fast;
    size_t depthsteps_max;
    double stop_step_incident;
    double stop_step_exiting;
    char *out_filename;
    char *exp_filename;
    char *bricks_filename;
    char *fit_vars;
    char *detector_out_filename;
    char *sample_filename;
    char *sample_out_filename;
    char **reaction_filenames;
    size_t n_reaction_filenames;
} cmdline_options;


#include "simulation.h"

void read_options(cmdline_options *cmd_opt, simulation *sim, int *argc, char ***argv);
cmdline_options *global_options_alloc();
void global_options_free(cmdline_options *cmd_opt);
const char *jabs_version();
void usage();

#endif // JABS_OPTIONS_H
