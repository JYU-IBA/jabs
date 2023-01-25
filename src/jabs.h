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
#include "geostragg.h"

double stop_sample(const sim_workspace *ws, const ion *incident, const sample *sample, gsto_stopping_type type, depth depth, double E);
depth next_crossing(const ion *incident, const sample *sample, const depth *d_from);
depth stop_step(const sim_workspace *ws, ion *incident, const sample *sample, struct depth depth, double step);
void exit_from_sample(ion *p, depth depth_start, const sim_workspace *ws, const sample *sample);
double stop_step_calculate(const sim_workspace *ws, const ion *ion);
int simulate(const ion *incident, depth depth_start, sim_workspace *ws, const sample *sample);
void simulate_reaction(sim_reaction *sim_r, const sim_workspace *ws, const sample *sample, const geostragg_vars *g, size_t i_depth, depth d_before, depth d_after, const ion *incident, double E_front, double S_front, double E_back, double S_back, double d_diff);
void simulate_init_reaction(sim_reaction *sim_r, const sample *sample, const geostragg_vars *g, double E_min, double E_max);
int assign_stopping(jibal_gsto *gsto, const simulation *sim);
int assign_stopping_Z2(jibal_gsto *gsto, const simulation *sim, int Z2); /* Assigns stopping and straggling (GSTO) for given Z2. Goes through all possible Z1s (beam and reaction products). */
int assign_stopping_Z1_Z2(jibal_gsto *gsto, int Z1, int Z2);
int print_spectra(const char *filename, const sim_workspace *ws, const gsl_histogram *exp);
int print_bricks(const char *filename, const sim_workspace *ws);
int simulate_with_ds(sim_workspace *ws);
double cross_section_concentration_product(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after, double S_front, double S_back);
double cross_section_concentration_product_fixed(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after, double S_front, double S_back);
double cross_section_concentration_product_adaptive(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after, double S_front, double S_back);
double cross_section_straggling(const sim_reaction *sim_r, gsl_integration_workspace *w, double accuracy, const prob_dist *pd, double E, double S);
double cross_section_straggling_fixed(const sim_reaction *sim_r, const prob_dist *pd, double E, double S);
double cross_section_straggling_adaptive(const sim_reaction *sim_r, gsl_integration_workspace *w, double accuracy, double E, double S);
fit_params *fit_params_all(fit_data *fit);
void fit_params_print(const fit_params *params, int active, const char *pattern); /* if active is TRUE print only active variables. pattern can be NULL to bypass matching. */
void fit_params_print_final(const fit_params *params);
size_t fit_params_enable(fit_params *params, const char *s, int enable); /* Enable/disable one or more variables matching pattern s. Returns number of matches. */
#endif // JABS_JABS_H
