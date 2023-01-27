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

typedef struct {
    depth d; /* Depth */
    double E; /* Energy */
    double S; /* Straggling */
} des; /* DES = Depth, Energy, Straggling */

typedef struct {
    des *t; /* array, n elements. E decreasing, d either increases (ion going deeper) or decreases. S can do whatever S does. */
    int depth_increases; /* TRUE if d increases (going deeper) */
    size_t n;
    size_t n_ranges;
    size_t *depth_interval_index; /* table, size same as number of sample ranges (as given in des_table_compute()). Array stores location i of t[i].d.x == sample->range[i_range].x  */
} des_table; /* Depth, energy, straggling table */

int des_table_realloc(des_table *dt, size_t n);
void des_table_free(des_table *dt);
des_table *des_table_compute(const ion *incident, depth depth_start, sim_workspace *ws, const sample *sample);
size_t des_table_size(const des_table *dt);
des *des_table_element(const des_table *dt, size_t i);
void des_table_rebuild_index(des_table *dt); /* called by des_table_compute() after setting values to table and before any other function can be used */
void des_table_print(FILE *f, const des_table *dt);
void des_table_set_ion_depth(const des_table *dt, ion *ion, depth d); /* Sets ion energy and straggling for some depth, based on des_table */
depth des_table_find_depth(const des_table *dt, size_t *i_des, depth depth_prev, ion *incident); /* Returns depth at given incident->E, or the next layer boundary, starting search from i_des in DES table. Updates i_des, incident->E and ->S. */
inline double des_table_min_energy(const des_table *dt) {return dt->t[dt->n - 1].E;}
double stop_sample(const sim_workspace *ws, const ion *incident, const sample *sample, gsto_stopping_type type, depth depth, double E);
depth next_crossing(const ion *incident, const sample *sample, const depth *d_from);
depth stop_step(const sim_workspace *ws, ion *incident, const sample *sample, struct depth depth, double step);
void exit_from_sample(ion *p, depth depth_start, const sim_workspace *ws, const sample *sample);
double stop_step_calculate(const sim_workspace *ws, const ion *ion);
double stop_step_calc_incident(const sim_workspace *ws, const ion *ion); /* New! */
int simulate(const ion *incident, depth depth_start, sim_workspace *ws, const sample *sample);
void simulate_reaction(sim_reaction *sim_r, const sim_workspace *ws, const sample *sample, const geostragg_vars *g, size_t i_depth, depth d_before, depth d_after, const ion *incident, double E_front, double S_front, double E_back, double S_back, double d_diff);
void simulate_reaction_new_routine(const ion *incident, depth depth_start, sim_workspace *ws, const sample *sample, const des_table *dt, const geostragg_vars *g, sim_reaction *sim_r);
void simulate_init_reaction(sim_reaction *sim_r, const sample *sample, const geostragg_vars *g, double E_min, double E_max);
int assign_stopping(jibal_gsto *gsto, const simulation *sim);
int assign_stopping_Z2(jibal_gsto *gsto, const simulation *sim, int Z2); /* Assigns stopping and straggling (GSTO) for given Z2. Goes through all possible Z1s (beam and reaction products). */
int assign_stopping_Z1_Z2(jibal_gsto *gsto, int Z1, int Z2);
int print_spectra(const char *filename, const sim_workspace *ws, const gsl_histogram *exp);
int print_bricks(const char *filename, const sim_workspace *ws);
int simulate_with_ds(sim_workspace *ws);
double cross_section_concentration_product_new(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after, double S_front, double S_back);
double cross_section_concentration_product(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after, double S_front, double S_back);
double cross_section_concentration_product_fixed(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after, double S_front, double S_back);
double cross_section_concentration_product_adaptive(const sim_workspace *ws, const sample *sample, const sim_reaction *sim_r, double E_front, double E_back, const depth *d_before, const depth *d_after, double S_front, double S_back);
double cross_section_straggling(const sim_reaction *sim_r, gsl_integration_workspace *w, double accuracy, const prob_dist *pd, double E, double S);
double cross_section_straggling_fixed(const sim_reaction *sim_r, const prob_dist *pd, double E, double S);
double cross_section_straggling_adaptive(const sim_reaction *sim_r, gsl_integration_workspace *w, double accuracy, double E, double S);

#endif // JABS_JABS_H
