/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2026 Jaakko Julin

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
#include "stop.h"
#include "des.h"

#ifdef __cplusplus
extern "C" {
#endif

int simulate(const ion *incident, depth depth_start, sim_workspace *ws, const sample *sample);
void simulate_reaction_ms(const sim_workspace *ws, const sim_reaction *sim_r, const ion *incident, double P_foil, brick *b); /* Calculates multiple scattering energy spread for brick b */
double ds_lateral_conc(const sim_workspace *ws, const sample *sample, size_t i_isotope, const ion *incident, double x); /* Concentration averaged over depth spread caused by lateral spread of ions in dual scattering */
int simulate_reaction(const ion *incident, depth depth_start, sim_workspace *ws, const sample *sample, const des_table *dt, const geostragg_vars *g, sim_reaction *sim_r);
int simulate_init_reaction(sim_reaction *sim_r, const sample *sample, const sim_calc_params *params, const geostragg_vars *g, double emin, double emin_incident, double emax_incident);
int assign_stopping(jibal_gsto *gsto, const simulation *sim);
int assign_stopping_Z2(jibal_gsto *gsto, const simulation *sim, int Z2); /* Assigns stopping and straggling (GSTO) for given Z2. Goes through all possible Z1s (beam and reaction products). */
int assign_stopping_Z1_Z2(jibal_gsto *gsto, int Z1, int Z2);
typedef struct ds_weight {
    double fluence; /* Fraction of incident ions scattered into this first scattering trajectory (polar angle bin and target isotope) */
    int skip; /* Trajectory is not calculated */
} ds_weight;

void ds_weights_skip_smallest(ds_weight *weights, size_t n, double fraction); /* Marks to be skipped the trajectories with the smallest weights, which together account for at most given fraction of total weight */
int simulate_with_ds(sim_workspace *ws);
double cross_section_concentration_product(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after, double S_front, double S_back);
double cross_section_concentration_product_adaptive(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after, double S_front, double S_back);
double cross_section_straggling(const sim_reaction *sim_r, gsl_integration_workspace *w, double accuracy, const prob_dist *pd, double E, double S, int emult);
double cross_section_straggling_fixed(const sim_reaction *sim_r, const prob_dist *pd, double E, double S, int emult);
double cross_section_straggling_adaptive(const sim_reaction *sim_r, gsl_integration_workspace *w, double accuracy, double E, double S, int emult);
double resonance_effect_mean_energy(const sim_reaction *sim_r, gsl_integration_workspace *w, double accuracy, const prob_dist *pd, double E, double S);
#ifdef __cplusplus
}
#endif
#endif // JABS_JABS_H
