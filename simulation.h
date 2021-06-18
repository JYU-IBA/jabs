/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_SIMULATION_H
#define JABS_SIMULATION_H

#include <gsl/gsl_histogram.h>
#include <jibal.h>
#include <time.h>

#include "ion.h"
#include "reaction.h"
#include "brick.h"
#include "detector.h"
#include "sample.h"

typedef struct {
    int ds;
    int ds_steps_azi;
    int ds_steps_polar;
    int n_ds;
    int cs_n_steps; /* Number of steps to take, when calculating cross section * concentration product */
    double cs_frac; /* Fractional step size 1.0/(cs_n_steps+1), calculated from cs_n_steps. */
    int cs_stragg_half_n;
    int cs_n_stragg_steps; /* calculated from cs_stragg_half_n */
    double emin;
    size_t depthsteps_max;
    int rk4;
    int nucl_stop_accurate;
    int mean_conc_and_energy;
    double stop_step_incident;
    double stop_step_exiting;
} sim_calc_params; /* All "calculation" parameters, i.e. not physical parameters */

typedef struct {
    reaction *reactions;
    size_t n_reactions;
    detector **det; /* Array of n_det detector pointers */
    size_t n_det;
    sample *sample;
    double p_sr;
    double sample_theta; /* Polar angle. Kind of. Zero is sample perpendicular to beam. */
    double sample_phi; /* Typically one uses a zero here, unless doing channeling stuff. Note that this is an azimuthal angle. */
    const jibal_isotope *beam_isotope;
    double beam_E;
    double beam_E_broad; /* Variance */
    double channeling_offset; /* a very ad-hoc channeling yield correction */
    double channeling_slope;
    sim_calc_params params;
} simulation;

typedef struct sim_reaction {
    const reaction *r;
    ion p; /* Reaction product */
    gsl_histogram *histo;
    brick *bricks;
    size_t n_bricks;
    size_t last_brick; /* inclusive, from 0 up to n_bricks-1 */
    int stop;
    double max_depth;
    size_t i_isotope; /* Number of isotope (r->target) in sample->isotopes */
    double theta;
    double K;
    double (*cross_section)(const struct sim_reaction *, double);
    double theta_cm; /* theta in CM system */
    double mass_ratio; /* m1/m2, these are variables to speed up cross section calculation */
    double E_cm_ratio;  /* m1/(m1+m2) */
    double cs_constant; /* Non-energy dependent Rutherford cross section terms for RBS or ERD */
    double r_VE_factor; /* Andersen correction factor r_VE = this / E_cm */
    double r_VE_factor2;
} sim_reaction; /* Workspace for a single reaction. Yes, the naming is confusing. */


typedef struct {
    double p_sr; /* With DS can be different from sim->p_sr, otherwise the same */
    const simulation *sim;
    const detector *det;
    const sample *sample; /* Note that simulate() can be passed a sample explicitly, but in most cases it should be this. Also this should be exactly the same as sim->sample. */
    size_t n_reactions;
    jibal_gsto *gsto;
    const jibal_config *jibal_config;
    gsto_stopping_type stopping_type;
    size_t n_channels; /* in histograms */
    gsl_histogram *histo_sum;
    ion ion;
    sim_reaction *reactions;
    const jibal_isotope *isotopes;
    sim_calc_params params;
} sim_workspace;


#include "sample.h"
simulation *sim_init();
void sim_free(simulation *sim);
sim_calc_params sim_calc_params_defaults(int ds, int fast);
int sim_reactions_add(simulation *sim, reaction_type type, jibal_cross_section_type cs, double theta); /* Add RBS or ERD reactions automagically */
int sim_sanity_check(const simulation *sim);
detector *sim_det(const simulation *sim, size_t i_det);
int sim_det_add(simulation *sim, detector *det);
int sim_det_set(simulation *sim, detector *det, size_t i_det); /* Will free existing detector (can be NULL too) */
sim_workspace *sim_workspace_init(const jibal *jibal, const simulation *sim, const detector *det);
void sim_workspace_free(sim_workspace *ws);
void sim_workspace_recalculate_n_channels(sim_workspace *ws, const simulation *sim);
void sim_workspace_calculate_sum_spectra(sim_workspace *ws);
void simulation_print(FILE *f, const simulation *sim);
void convolute_bricks(sim_workspace *ws);
void sim_reaction_recalculate_internal_variables(sim_reaction *sim_r);
void sim_reaction_reset_bricks(sim_reaction *sim_r);
double sim_reaction_cross_section_rutherford(const sim_reaction *sim_r, double E);
double sim_reaction_cross_section_tabulated(const sim_reaction *sim_r, double E);
double sim_reaction_andersen(const sim_reaction *sim_r, double E_cm);
#endif // JABS_SIMULATION_H
