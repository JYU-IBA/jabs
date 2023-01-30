/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include "defaults.h"
#include "message.h"
#include "sim_calc_params.h"

sim_calc_params *sim_calc_params_defaults(sim_calc_params *p) {
    if(!p) {
        p = malloc(sizeof(sim_calc_params));
        p->cs_stragg_pd = NULL; /* Will be allocated by sim_calc_params_update() */
    }
    p->incident_stop_params.step = INCIDENT_STOP_STEP_DEFAULT;
    p->incident_stop_params.sigmas = INCIDENT_STOP_STEP_SIGMAS_DEFAULT;
    p->incident_stop_params.min = INCIDENT_STOP_STEP_MIN_DEFAULT;
    p->incident_stop_params.max = INCIDENT_STOP_STEP_MAX_DEFAULT;
    p->exiting_stop_params.step = EXITING_STOP_STEP_DEFAULT;
    p->exiting_stop_params.sigmas = EXITING_STOP_STEP_SIGMAS_DEFAULT;
    p->exiting_stop_params.min = EXITING_STOP_STEP_MIN_DEFAULT;
    p->exiting_stop_params.max = EXITING_STOP_STEP_MAX_DEFAULT;
    p->brick_width_sigmas = BRICK_WIDTH_SIGMAS_DEFAULT;
    p->n_bricks_max = 0; /* automatic */
    p->geostragg = FALSE;
    p->beta_manual  = FALSE;
    p->ds = FALSE;
    p->ds_steps_azi = 0;
    p->ds_steps_polar = 0;
    p->rk4 = TRUE;
    p->nuclear_stopping_accurate = TRUE;
    p->mean_conc_and_energy = FALSE;
    p->cs_n_stragg_steps = CS_STRAGG_STEPS;
    p->rough_layer_multiplier = 1.0;
    p->sigmas_cutoff = SIGMAS_CUTOFF;
    p->gaussian_accurate = FALSE;
    p->int_cs_max_intervals = CS_CONC_MAX_INTEGRATION_INTERVALS;
    p->int_cs_accuracy = CS_CONC_INTEGRATION_ACCURACY;
    p->int_cs_stragg_max_intervals = CS_STRAGG_MAX_INTEGRATION_INTERVALS;
    p->int_cs_stragg_accuracy = CS_STRAGG_INTEGRATION_ACCURACY;
    p->cs_adaptive = FALSE;
    p->cs_energy_step_max = CS_ENERGY_STEP_MAX_DEFAULT;
    p->cs_depth_step_max = CS_DEPTH_STEP_MAX_DEFAULT;
    p->cs_stragg_step_sigmas = CS_STRAGG_STEP_FUDGE_FACTOR_DEFAULT;
    p->ds_incident_stop_step_factor = DUAL_SCATTER_INCIDENT_STOP_STEP_FACTOR_DEFAULT;
    p->reaction_file_angle_tolerance = REACTION_FILE_ANGLE_TOLERANCE_DEFAULT;
#ifdef DEBUG
    fprintf(stderr, "New calc params created.\n");
#endif
    return p;
}

sim_calc_params *sim_calc_params_defaults_fast(sim_calc_params *p) {
    sim_calc_params_defaults(p);
    p->rk4 = FALSE;
    p->nuclear_stopping_accurate = FALSE;
    p->mean_conc_and_energy = TRUE;
    p->cs_n_stragg_steps = 0; /* Not used if mean_conc_and_energy == TRUE */
    p->geostragg = FALSE;
    p->rough_layer_multiplier = 0.5;
    p->sigmas_cutoff = SIGMAS_FAST_CUTOFF;
    p->incident_stop_params.min *= 2.0;
    p->incident_stop_params.sigmas *= 1.5;
    p->exiting_stop_params.min *= 1.5;
    p->exiting_stop_params.sigmas *= 1.5;
    p->exiting_stop_params.max *= 1.5;
    return p;
}

sim_calc_params *sim_calc_params_defaults_accurate(sim_calc_params *p) {
    sim_calc_params_defaults(p);
    p->cs_n_stragg_steps = 0; /* Automatic (adaptive) */
    p->sigmas_cutoff += 1.0;
    p->gaussian_accurate = TRUE;
    p->cs_adaptive = TRUE;
    p->incident_stop_params.min *= 0.5;
    p->exiting_stop_params.min *= 0.5;
    p->exiting_stop_params.sigmas *= 0.5;
    return p;
}

sim_calc_params *sim_calc_params_defaults_brisk(sim_calc_params *p) {
    sim_calc_params_defaults(p);
    p->cs_n_stragg_steps -= 2;
    p->sigmas_cutoff -= 1.0;
    p->incident_stop_params.min *= 2.0;
    p->exiting_stop_params.sigmas *= 1.5;
    p->cs_energy_step_max *= 1.5;
    p->cs_depth_step_max *= 1.5;
    p->cs_stragg_step_sigmas = 1.25;
    return p;
}

sim_calc_params *sim_calc_params_defaults_improved(sim_calc_params *p) {
    sim_calc_params_defaults(p);
    p->cs_n_stragg_steps += 4;
    p->sigmas_cutoff += 0.5;
    p->gaussian_accurate = TRUE;
    p->exiting_stop_params.min *= 0.75;
    p->exiting_stop_params.sigmas *= 0.75;
    p->cs_energy_step_max *= 0.75;
    p->cs_depth_step_max *= 0.75;
    p->cs_stragg_step_sigmas = 0.75;
    return p;
}

void sim_calc_params_free(sim_calc_params *p) {
    if(!p)
        return;
    prob_dist_free(p->cs_stragg_pd);
    free(p);
}

void sim_calc_params_copy(const sim_calc_params *p_src, sim_calc_params *p_dst) {
    *p_dst = *p_src;
    p_dst->cs_stragg_pd = NULL;
}

void sim_calc_params_update(sim_calc_params *p) {
    prob_dist_free(p->cs_stragg_pd);
    if(p->mean_conc_and_energy) {
        p->cs_n_stragg_steps = 0;
    }
    p->cs_stragg_pd = prob_dist_gaussian(p->cs_n_stragg_steps);

    if(p->ds && p->ds_steps_azi == 0 && p->ds_steps_polar == 0) { /* DS defaults are applied if nothing else is specified */
        p->ds_steps_azi = DUAL_SCATTER_AZI_STEPS;
        p->ds_steps_polar = DUAL_SCATTER_POLAR_STEPS;
    }
}

void sim_calc_params_ds(sim_calc_params *p, int ds) {
    if(ds) {
        p->ds = TRUE;
    }
}

void sim_calc_params_print(const sim_calc_params *params) {
    if(!params)
        return;
    jabs_message(MSG_INFO, stderr, "step for incident ions = %.3lf keV (0 = auto)\n", params->incident_stop_params.step / C_KEV);
    if(params->incident_stop_params.step == 0.0) {
        jabs_message(MSG_INFO, stderr, "step for incident ions = %g times straggling sigma\n", params->incident_stop_params.sigmas);
        jabs_message(MSG_INFO, stderr, "minimum step for incident ions = %.3lf keV\n", params->incident_stop_params.min / C_KEV);
        jabs_message(MSG_INFO, stderr, "maximum step for incident ions = %.3lf keV\n", params->incident_stop_params.max / C_KEV);
    }
    jabs_message(MSG_INFO, stderr, "step for exiting ions = %.3lf keV (0 = auto)\n", params->exiting_stop_params.step / C_KEV);
    if(params->exiting_stop_params.step == 0.0) {
        jabs_message(MSG_INFO, stderr, "step for exiting ions = %g times straggling sigma\n", params->exiting_stop_params.sigmas);
        jabs_message(MSG_INFO, stderr, "minimum step for exiting ions = %.3lf keV\n", params->exiting_stop_params.min / C_KEV);
        jabs_message(MSG_INFO, stderr, "maximum step for exiting ions = %.3lf keV\n", params->exiting_stop_params.max / C_KEV);
    }
    jabs_message(MSG_INFO, stderr, "stopping RK4 = %s\n", params->rk4?"true":"false");
    jabs_message(MSG_INFO, stderr, "accurate nuclear stopping = %s\n", params->nuclear_stopping_accurate?"true":"false");
    jabs_message(MSG_INFO, stderr, "maximum number of bricks = %zu\n", params->n_bricks_max);
    jabs_message(MSG_INFO, stderr, "geometric broadening = %s\n", params->geostragg?"true":"false");
    jabs_message(MSG_INFO, stderr, "brick width = %g times detector and straggling sum sigma\n", params->brick_width_sigmas);
    jabs_message(MSG_INFO, stderr, "cross section of brick determined using mean concentration and energy = %s\n", params->mean_conc_and_energy?"true":"false");
    if(!params->mean_conc_and_energy) {
        if(params->cs_adaptive) {
            jabs_message(MSG_INFO, stderr, "cross section integration accuracy = %g\n", params->int_cs_accuracy);
        } else {
            jabs_message(MSG_INFO, stderr, "cross section evaluation step = %g times straggling sigma\n", params->cs_stragg_step_sigmas);
            jabs_message(MSG_INFO, stderr, "maximum energy step for cross section evaluation = %g keV\n", params->cs_energy_step_max / C_KEV);
            jabs_message(MSG_INFO, stderr, "maximum depth step for cross section evaluation = %g tfu\n", params->cs_depth_step_max / C_TFU);
        }
        if(params->cs_n_stragg_steps == 0) {
            jabs_message(MSG_INFO, stderr, "straggling weighting integration accuracy = %g\n", params->int_cs_stragg_accuracy);
        } else {
            jabs_message(MSG_INFO, stderr, "straggling substeps = %zu\n", params->cs_n_stragg_steps);
        }
    }
    jabs_message(MSG_INFO, stderr, "reaction file (R33) angle tolerance = %g deg\n", params->reaction_file_angle_tolerance / C_DEG);
}
