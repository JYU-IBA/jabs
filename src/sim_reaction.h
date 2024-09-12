/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2024 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_SIM_REACTION_H
#define JABS_SIM_REACTION_H

#include "reaction.h"
#include "ion.h"
#include "brick.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct sim_reaction {
    const reaction *r;
    ion p; /* Reaction product */
    jabs_histogram *histo;
    brick *bricks;
    size_t n_bricks;
    size_t last_brick; /* inclusive, from 0 up to n_bricks-1 */
    int stop;
    size_t n_convolution_calls; /* Set to zero on init, but not reset after that. Incremented by sim_workspace_histograms_calculate(). Nonzero means a spectrum probably exists. */
    double max_depth;
    size_t i_isotope; /* Number of isotope (r->target) in sample->isotopes */
    double theta; /* theta used in simulations (will be overwritten by simulate() when necessary) */
    double K;
    double (*cross_section)(const struct sim_reaction *, double);
    double theta_cm; /* theta in CM system */
    double mass_ratio; /* m1/m2, these are variables to speed up cross section calculation */
    double E_cm_ratio;  /* m1/(m1+m2) */
    double cs_constant; /* Non-energy dependent Rutherford cross section terms for RBS or ERD */
    double r_VE_factor; /* Andersen correction factor r_VE = this / E_cm */
    double r_VE_factor2;
    double lecuyer_factor; /* F_LEcyuer = 1 - this / E_cm */
    double emin_incident; /* Set on sim_reaction_recalculate_internal_variables(), should reflect lowest usable TODO: *incident*? energy  by "ws", reaction and stopping */
    double emax_incident; /* Set on sim_reaction_recalculate_internal_variables(), should reflect highest usable TODO: *incident*? energy limited by "ws", reaction, stopping and ion (incident E_0) */
    double emin_product;
    double emax_product;
    struct reaction_point *cs_table; /* precalculated screening corrections from emin to emax, n_cs_table points */
    size_t n_cs_table;
    double cs_estep; /* calculated by sim_reaction_recalculate_screening_table(), step size based on n_cs_table */
} sim_reaction; /* Workspace for a single reaction. Yes, the naming is confusing. */

sim_reaction *sim_reaction_init(const sample *sample, const detector *det, const reaction *r, size_t n_channels, size_t n_bricks);
void sim_reaction_free(sim_reaction *sim_r);
int sim_reaction_recalculate_internal_variables(sim_reaction *sim_r, const sim_calc_params *params, double theta, double emin, double emin_incident, double emax_incident);
int sim_reaction_recalculate_screening_table(sim_reaction *sim_r);
void sim_reaction_reset_screening_table(sim_reaction *sim_r);
double sim_reaction_evaluate_screening_table(const sim_reaction *sim_r, double E);
void sim_reaction_reset_bricks(sim_reaction *sim_r);
void sim_reaction_set_cross_section_by_type(sim_reaction *sim_r);
double sim_reaction_cross_section_rutherford(const sim_reaction *sim_r, double E);
double sim_reaction_cross_section_tabulated(const sim_reaction *sim_r, double E);
#ifdef JABS_PLUGINS
double sim_reaction_cross_section_plugin(const sim_reaction *sim_r, double E);
#endif
double sim_reaction_andersen(const sim_reaction *sim_r, double E_cm);
double sim_reaction_lecuyer(const sim_reaction *sim_r, double E_cm);
void sim_reaction_product_energy_and_straggling(sim_reaction *r, const ion *incident);
void sim_reaction_print_bricks(FILE *f, const sim_reaction *r, double psr);
#ifdef __cplusplus
}
#endif
#endif // JABS_SIM_REACTION_H
