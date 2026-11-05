/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2026 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_SIMULATION_WORKSPACE_H
#define JABS_SIMULATION_WORKSPACE_H
#include <gsl/gsl_integration.h>
#include "simulation.h"
#include "sim_reaction.h"
#include "spectrum.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct ds_trajectory { /* Dual scattering trajectory information, used for path length criterion of trajectories with a small scattering angle (Barradas, NIMB 225 (2004) 318) */
    int active; /* TRUE when simulating dual scattering (second scattering) and criterion should be checked */
    int first_small; /* First scattering angle is small */
    double x1; /* Depth of first scattering */
    double l1; /* Path length of incident ion before first scattering */
    double inverse_cosine_in; /* Absolute value of inverse cosine of incident beam (in sample) */
    double angle_min; /* Scattering angles below this are small */
    double path_cutoff; /* Path length ratio required */
    int lateral; /* TRUE when lateral spread (depth spread) of ions between first and second scattering is taken into account */
} ds_trajectory;

typedef struct sim_workspace {
    double fluence; /* With DS can be different from sim->fluence, otherwise the same */
    const simulation *sim;
    const detector *det;
    const sample *sample; /* Note that simulate() can be passed a sample explicitly, but in most cases it should be this. Also this should be exactly the same as sim->sample. */
    size_t n_reactions;
    const jibal_gsto *gsto;
    size_t n_channels; /* in histograms */
    jabs_histogram *histo_sum;
    jabs_histogram *histo_ds; /* Dual scattering contribution to histo_sum (sum of reaction histo_ds), NULL if dual scattering is not calculated */
    int ds_running; /* TRUE while dual scattering trajectories are calculated, spectra are then convoluted to histo_ds of reactions */
    ion ion;
    sim_reaction **reactions; /* table of reaction pointers, size n_reactions */
    const jibal_isotope *isotopes;
    sim_calc_params *params;
    jabs_stop stop; /* Stopping calculation parameters and data, set on sim_workspace_init() */
    jabs_stop stragg; /* Straggling calculation parameters and data */
    ms_model *ms; /* Multiple scattering model, NULL if not needed */
    int ms_spread; /* Energy spread due to multiple scattering is calculated (TRUE/FALSE) */
    ds_trajectory ds_traj;
    double ds_beam_loss_angle; /* If larger than zero, incident beam is attenuated by scattering to larger angles than this (dual scattering) */
    double emin;
    gsl_integration_workspace *w_int_cs; /* Integration workspace for conc * cross section product */
    gsl_integration_workspace *w_int_cs_stragg;
    size_t n_bricks; /* same as r->n_bricks in each reaction */
} sim_workspace;


sim_workspace *sim_workspace_init(const jibal *jibal, const simulation *sim, const detector *det);
void sim_workspace_init_reactions(sim_workspace *ws); /* used by sim_workspace_init(), ws->sim and ws->n_bricks should be set before calling */
void sim_workspace_calculate_number_of_bricks(sim_workspace *ws);
void sim_workspace_free(sim_workspace *ws);
void sim_workspace_recalculate_n_channels(sim_workspace *ws, const simulation *sim);
void sim_workspace_calculate_sum_spectra(sim_workspace *ws);

void sim_workspace_histograms_reset(sim_workspace *ws);
size_t sim_workspace_histograms_calculate(sim_workspace *ws);
void sim_workspace_histograms_scale(sim_workspace *ws, double scale);
int sim_workspace_print_spectra(const result_spectra *spectra, const char *filename, int cl);
int sim_workspace_print_bricks(const sim_workspace *ws, const char *filename);
#ifdef __cplusplus
}
#endif
#endif //JABS_SIMULATION_WORKSPACE_H
