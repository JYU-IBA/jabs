/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_JABS_H
#define JABS_JABS_H

#include "simulation.h"
#include "ion.h"
#include "fit.h"
#include "options.h"
#include "sample.h"
#include "reaction.h"

double stop_sample(const sim_workspace *ws, const ion *incident, const sample *sample, gsto_stopping_type type, depth depth, double E);
depth next_crossing(const ion *incident, const sample *sample, const depth *d_from);
depth stop_step(const sim_workspace *ws, ion *incident, const sample *sample, struct depth depth, double step);
void post_scatter_exit(ion *p, depth depth_start, const sim_workspace *ws, const sample *sample);
void foil_traverse(ion *p, const sample *foil, sim_workspace *ws);
int simulate(const ion *incident, depth depth_start, sim_workspace *ws, const sample *sample);
int assign_stopping(jibal_gsto *gsto, const simulation *sim);
int assign_stopping_Z2(jibal_gsto *gsto, const simulation *sim, int Z2); /* Assigns stopping and straggling (GSTO) for given Z2. Goes through all possible Z1s (beam and reaction products). */
int assign_stopping_Z1_Z2(jibal_gsto *gsto, int Z1, int Z2);
int print_spectra(const char *filename, const sim_workspace *ws, const gsl_histogram *exp);
int fit_params_add(simulation *sim, const sample_model *sm, fit_params *params, const char *fit_vars); /* fit_vars is a comma separated list of variables to fit */
void fit_params_clear(fit_params *params);
int print_bricks(const char *filename, const sim_workspace *ws);
int simulate_with_ds(sim_workspace *ws);
void ds(sim_workspace *ws); /* TODO: the DS routine is more pseudocode at this stage... */
double cross_section_concentration_product(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after, double S_front, double S_back);
#endif // JABS_JABS_H
