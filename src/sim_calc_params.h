/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2026 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_SIM_CALC_PARAMS_H
#define JABS_SIM_CALC_PARAMS_H
#include "stop.h"
#include "prob_dist.h"
#include "message.h"
#include "ms.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum ds_reaction_spectra_type {
    DS_REACTION_SPECTRA_INCLUDED = 0, /* Dual scattering is included in reaction spectra (they sum up to the simulated spectrum) */
    DS_REACTION_SPECTRA_SEPARATE = 1, /* Reaction spectra are single scattering only, dual scattering of each reaction is output separately */
    DS_REACTION_SPECTRA_NONE = 2 /* Reaction spectra are single scattering only, dual scattering is only in the sum spectrum */
} ds_reaction_spectra_type;

static const jibal_option ds_reaction_spectra_option[] = {
        {"included", DS_REACTION_SPECTRA_INCLUDED},
        {"separate", DS_REACTION_SPECTRA_SEPARATE},
        {JIBAL_OPTION_STR_NONE, DS_REACTION_SPECTRA_NONE},
        {NULL, 0}
};

typedef struct sim_calc_params {
    int ds; /* Dual scattering true/false */
    int ds_steps_azi;
    int ds_steps_polar;
    int ds_steps_polar_substeps;
    int ds_small_angles; /* Include trajectories with one small scattering angle (TRUE/FALSE) */
    double ds_small_angle_min; /* Smallest scattering angle considered when ds_small_angles is TRUE */
    int ds_small_angle_steps; /* Polar steps between ds_small_angle_min and DUAL_SCATTER_ANGLE_MIN */
    double ds_path_cutoff; /* Path length ratio criterion for trajectories with a small scattering angle */
    int ds_beam_loss; /* Incident beam is attenuated by dual scattering (TRUE/FALSE) */
    int ds_lateral_spread; /* Lateral spread due to multiple scattering is taken into account between first and second scattering (TRUE/FALSE) */
    double ds_skip_fraction; /* Dual scattering trajectories with smallest weights that account for this fraction of total weight (at each depth) are not calculated */
    double ds_grazing_angle; /* Dual scattering trajectories (after first scattering) closer than this angle to the sample surface are not calculated */
    ds_reaction_spectra_type ds_reaction_spectra; /* How dual scattering is shown in reaction spectra (output only). Sum of dual scattering is always output separately. */
    int ds_accurate; /* Dual scattering trajectories are calculated with the same parameters as single scattering (TRUE) or with parameters adjusted for speed (FALSE, see sim_calc_params_ds_trajectories()) */
    prob_dist *cs_stragg_pd;
    size_t cs_n_stragg_steps; /* Number of steps to take, when calculating straggling weighted cross sections. If zero, adaptive integration is used. */
    size_t n_bricks_max;
    int rk4; /* Use fourth order Runge-Kutta for energy loss calculation (differential equation with dE/dx). When false, a first-order method is used. */
    int nuclear_stopping_accurate; /* Use accurate nuclear stopping equation true/false. When false a faster (poorly approximating) equation is used below the nuclear stopping maximum. */
    int mean_conc_and_energy; /* Calculation of cross-section concentration product is simplified by calculating cross section at mean energy of a depth step and concentration at mid-bin (only relevant for samples with concentration gradients) */
    int geostragg; /* Geometric straggling true/false */
    ms_model_type ms; /* Small angle multiple scattering model, MS_MODEL_NONE if disabled */
    int beta_manual; /* Don't calculate exit angle based on detector geometry, use something given by user, true/false */
    int gaussian_accurate; /* If this is FALSE an approximative gaussian CDF is used in convolution of spectra, otherwise function from GSL is used. */
    jabs_stop_step_params incident_stop_params;
    jabs_stop_step_params exiting_stop_params;
    double ds_incident_stop_step_factor;
    double brick_width_sigmas;
    double rough_layer_multiplier; /* Multiply given (or default) number of subspectra when calculating rough layers. */
    double sigmas_cutoff; /* Number of (+-) sigmas to consider when turning bricks to spectra */
    size_t int_cs_max_intervals;
    double int_cs_accuracy; /* Accuracy of conc*cross section integration (not always relevant) */
    size_t int_cs_stragg_max_intervals;
    double int_cs_stragg_accuracy; /* Accuracy of conc * straggling (gaussian) integration (not always relevant) */
    int cs_adaptive;
    double cs_energy_step_max; /* Largest energy step (incident beam, mean) between cross section * concentration evaluation. Weight by straggling happens inside this step (finer stepping).*/
    double cs_depth_step_max /* Largest depth step for section * concentration evaluation. Also see cs_energy_step_max. */ ;
    double cs_stragg_step_sigmas; /*  Cross section * concentration evaluation straggling step multiplier. Step is standard deviation times this. Maximum step may also be limited by cs_depth_step_max, cs_energy_step_max  */
    double reaction_file_angle_tolerance;
    int bricks_skip_zero_conc_ranges; /* TRUE/FALSE, determines if we should calculate so-called "empty" bricks or skip over them.*/
    int screening_tables;
    int bricks_save; /* Saving of brick data to text files (for debugging purposes etc) */
    int mean_energy_after_reaction_accurate; /* Calculate mean energy after reaction with for all reactions, not for REACTION_FILE type reactions. */
} sim_calc_params; /* All "calculation" parameters, i.e. not physical parameters */

sim_calc_params *sim_calc_params_defaults(sim_calc_params *p); /* if p is NULL, allocates params */
sim_calc_params *sim_calc_params_defaults_fast(sim_calc_params *p); /* Sets parameters to defaults and then makes them faster */
sim_calc_params *sim_calc_params_defaults_accurate(sim_calc_params *p); /* Sets parameters to defaults and then makes them a lot slower */
sim_calc_params *sim_calc_params_defaults_brisk(sim_calc_params *p); /* Sets parameters to defaults and then makes them slightly faster */
sim_calc_params *sim_calc_params_defaults_improved(sim_calc_params *p); /* Sets parameters to defaults and then makes them slightly faster */
void sim_calc_params_free(sim_calc_params *p);
void sim_calc_params_copy(const sim_calc_params *p_src, sim_calc_params *p_dst);
void sim_calc_params_update(sim_calc_params *p); /* Computes variables that can be computed from other variables */
void sim_calc_params_ds(sim_calc_params *p, int ds); /* if ds is TRUE, set DS parameters, otherwise no action is taken */
void sim_calc_params_ds_trajectories(sim_calc_params *p); /* Adjusts (current) parameters for calculation of dual scattering trajectories, unless p->ds_accurate is TRUE. DS settings are not changed. */
void sim_calc_params_print(const sim_calc_params *params, jabs_msg_level msg_level);
#ifdef __cplusplus
}
#endif
#endif //JABS_SIM_CALC_PARAMS_H
