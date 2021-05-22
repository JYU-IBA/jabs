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

double stop_sample(sim_workspace *ws, const ion *incident, const sample *sample, gsto_stopping_type type, struct depth depth, double E);
depth next_crossing(const ion *incident, const sample *sample, const depth *d_from);
depth stop_step(sim_workspace *ws, ion *incident, const sample *sample, struct depth depth, double step);
void simulate(const ion *incident, depth depth_start, sim_workspace *ws, const sample *sample);
reaction **make_reactions(const sample *sample, const simulation *sim, jibal_cross_section_type cs_rbs, jibal_cross_section_type cs_erd);/* Note that sim->ion needs to be set and geometry has to be correct */
int process_reaction_files(const jibal_isotope *jibal_isotopes, reaction **reactions, char * const *reaction_filenames, size_t n_reaction_filenames);
int assign_stopping(jibal_gsto *gsto, const simulation *sim, const sample *sample, reaction * const *reactions);
int print_spectra(const char *filename, const sim_workspace *ws, const gsl_histogram *exp);
void add_fit_params(simulation *sim, const sample_model *sm, fit_params *params, const char *fit_vars); /* fit_vars is a comma separated list of variables to fit */
void output_bricks(const char *filename, const sim_workspace *ws);
void simulate_with_ds(sim_workspace *ws);
void ds(sim_workspace *ws); /* TODO: the DS routine is more pseudocode at this stage... */
double cross_section_concentration_product(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after, double S_front, double S_back);
#endif // JABS_JABS_H
