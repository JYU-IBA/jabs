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
#include <gsl/gsl_integration.h>
#include <jibal.h>
#include <time.h>

#include "ion.h"
#include "reaction.h"
#include "brick.h"
#include "detector.h"
#include "sample.h"
#include "prob_dist.h"


typedef struct {
    int ds; /* Dual scattering true/false */
    int ds_steps_azi;
    int ds_steps_polar;
    int n_ds;
    size_t cs_n_steps; /* Number of steps to take, when calculating cross section * concentration product */
    double cs_frac; /* Fractional step size 1.0/(cs_n_steps+1), calculated from cs_n_steps. */
    prob_dist *cs_stragg_pd;
    size_t cs_n_stragg_steps; /* Number of steps to take, when calculating straggling weighted cross sections (note that these are substeps of cs_n_steps) */
    size_t depthsteps_max;
    int rk4; /* Use fourth order Runge-Kutta for energy loss calculation (differential equation with dE/dx). When false, a first-order method is used. */
    int nucl_stop_accurate; /* Use accurate nuclear stopping equation true/false. When false a faster (poorly approximating) equation is used below the nuclear stopping maximum. */
    int mean_conc_and_energy; /* Calculation of cross-section concentration product is simplified by calculating cross section at mean energy of a depth step and concentration at mid-bin (only relevant for samples with concentration gradients) */
    int geostragg; /* Geometric straggling true/false */
    int beta_manual; /* Don't calculate exit angle based on detector geometry, use something given by user, true/false */
    int gaussian_accurate; /* If this is FALSE an approximative gaussian CDF is used in convolution of spectra, otherwise function from GSL is used. */
    double stop_step_incident;
    double stop_step_exiting;
    double stop_step_fudge_factor;
    double stop_step_min; /* TODO: automatic */
    double stop_step_add; /* This is added to stop step */
    double rough_layer_multiplier; /* Multiply given (or default) number of subspectra when calculating rough layers. */
    double sigmas_cutoff; /* Number of (+-) sigmas to consider when turning bricks to spectra */
    size_t int_cs_max_intervals;
    double int_cs_accuracy; /* Accuracy of conc*cross section integration (not always relevant) */
    size_t int_cs_stragg_max_intervals;
    double int_cs_stragg_accuracy; /* Accuracy of conc * straggling (gaussian) integration (not always relevant) */
} sim_calc_params; /* All "calculation" parameters, i.e. not physical parameters */

typedef struct {
    reaction **reactions;
    size_t n_reactions;
    detector **det; /* Array of n_det detector pointers */
    size_t n_det;
    sample *sample;
    double fluence;
    double sample_theta; /* Polar angle. Kind of. Zero is sample perpendicular to beam. */
    double sample_phi; /* Typically one uses a zero here, unless doing channeling stuff. Note that this is an azimuthal angle. */
    const jibal_isotope *beam_isotope;
    aperture *beam_aperture;
    double beam_E;
    double beam_E_broad; /* Variance */
    double emin;
    double channeling_offset; /* a very ad-hoc channeling yield correction */
    double channeling_slope;
    sim_calc_params *params;
    int erd; /* Add ERD reactions */
    int rbs; /* Add RBS reactions */
    jibal_cross_section_type cs_rbs;
    jibal_cross_section_type cs_erd;
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
    double fluence; /* With DS can be different from sim->fluence, otherwise the same */
    const simulation *sim;
    const detector *det;
    const sample *sample; /* Note that simulate() can be passed a sample explicitly, but in most cases it should be this. Also this should be exactly the same as sim->sample. */
    size_t n_reactions;
    jibal_gsto *gsto;
    gsto_stopping_type stopping_type;
    size_t n_channels; /* in histograms */
    gsl_histogram *histo_sum;
    ion ion;
    sim_reaction *reactions;
    const jibal_isotope *isotopes;
    sim_calc_params *params;
    double emin;
    gsl_integration_workspace *w_int_cs; /* Integration workspace for conc * cross section product */
    gsl_integration_workspace *w_int_cs_stragg;
} sim_workspace;


#include "sample.h"
simulation *sim_init(jibal *jibal);
void sim_free(simulation *sim);
sim_calc_params *sim_calc_params_defaults(sim_calc_params *p); /* if p is NULL, allocates params */
sim_calc_params *sim_calc_params_defaults_fast(sim_calc_params *p); /* Sets parameters to defaults and then makes them faster */
sim_calc_params *sim_calc_params_defaults_accurate(sim_calc_params *p); /* Sets parameters to defaults and then makes them a lot slower */
sim_calc_params *sim_calc_params_defaults_brisk(sim_calc_params *p); /* Sets parameters to defaults and then makes them slightly faster */
sim_calc_params *sim_calc_params_defaults_improved(sim_calc_params *p); /* Sets parameters to defaults and then makes them slightly faster */
void sim_calc_params_free(sim_calc_params *p);
void sim_calc_params_copy(const sim_calc_params *p_src, sim_calc_params *p_dst);
void sim_calc_params_update(sim_calc_params *p); /* Computes variables that can be computed from other variables */
void sim_calc_params_ds(sim_calc_params *p, int ds); /* if ds is TRUE, set DS parameters, otherwise no action is taken */
void sim_calc_params_print(const sim_calc_params *params);
jibal_cross_section_type sim_cs(const simulation *sim, reaction_type type);
int sim_reactions_add_reaction(simulation *sim, reaction *r);
int sim_reactions_remove_reaction(simulation *sim, size_t i);
int sim_reactions_add_auto(simulation *sim, const sample_model *sm, reaction_type type, jibal_cross_section_type cs); /* Add RBS or ERD reactions automagically */
int sim_reactions_add_r33(simulation *sim, const jibal_isotope *jibal_isotopes, const char *filename);
void sim_reactions_free(simulation *sim); /* Free reactions and reset the number of reactions to zero */
int sim_sanity_check(const simulation *sim);
detector *sim_det(const simulation *sim, size_t i_det);
detector *sim_det_from_string(const simulation *sim, const char *s);
int sim_det_add(simulation *sim, detector *det);
int sim_det_set(simulation *sim, detector *det, size_t i_det); /* Will free existing detector (can be NULL too) */
sim_workspace *sim_workspace_init(const jibal *jibal, const simulation *sim, const detector *det);
void sim_workspace_free(sim_workspace *ws);
void sim_workspace_recalculate_n_channels(sim_workspace *ws, const simulation *sim);
void sim_workspace_calculate_sum_spectra(sim_workspace *ws);
void sim_print(const simulation *sim);
void sim_workspace_histograms_reset(sim_workspace *ws);
void sim_workspace_histograms_calculate(sim_workspace *ws);
void sim_workspace_histograms_scale(sim_workspace *ws, double scale);
void sim_reaction_recalculate_internal_variables(sim_reaction *sim_r);
void sim_reaction_reset_bricks(sim_reaction *sim_r);
double sim_reaction_cross_section_rutherford(const sim_reaction *sim_r, double E);
double sim_reaction_cross_section_tabulated(const sim_reaction *sim_r, double E);
double sim_reaction_andersen(const sim_reaction *sim_r, double E_cm);
void sim_sort_reactions(const simulation *sim);
void sim_reaction_product_energy_and_straggling(sim_reaction *r, const ion *incident);
double sim_alpha_angle(const simulation *sim);
double sim_exit_angle(const simulation *sim, const detector *det);
int sim_do_we_need_erd(const simulation *sim);
#endif // JABS_SIMULATION_H
