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

double stop_sample(sim_workspace *ws, const ion *incident, const sample *sample, gsto_stopping_type type, const depth depth, double E);
depth stop_step(sim_workspace *ws, ion *incident, const sample *sample, const depth depth, double step);
void simulate(const ion *incident, double x_0, sim_workspace *ws, const sample *sample);
reaction make_reaction(const sample *sample, const simulation *sim, size_t i_isotope, reaction_type type);
reaction *make_reactions(const sample *sample, const simulation *sim, int rbs, int erd);/* Note that sim->ion needs to be set and geometry has to be correct */
int assign_stopping(jibal_gsto *gsto, const simulation *sim, const sample *sample, const reaction *reactions);
void print_spectra(FILE *f, const global_options *global, const sim_workspace *ws, const gsl_histogram *exp);
void add_fit_params(global_options *global, simulation *sim, jibal_layer **layers, size_t n_layers, fit_params *params);
void output_bricks(const char *filename, const sim_workspace *ws);
void no_ds(sim_workspace *ws, const sample *sample);
void ds(sim_workspace *ws, const sample *sample); /* TODO: the DS routine is more pseudocode at this stage... */

#endif // JABS_JABS_H
