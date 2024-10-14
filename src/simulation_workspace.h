/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2024 Jaakko Julin

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


typedef struct sim_workspace {
    double fluence; /* With DS can be different from sim->fluence, otherwise the same */
    const simulation *sim;
    const detector *det;
    const sample *sample; /* Note that simulate() can be passed a sample explicitly, but in most cases it should be this. Also this should be exactly the same as sim->sample. */
    size_t n_reactions;
    const jibal_gsto *gsto;
    size_t n_channels; /* in histograms */
    jabs_histogram *histo_sum;
    ion ion;
    sim_reaction **reactions; /* table of reaction pointers, size n_reactions */
    const jibal_isotope *isotopes;
    sim_calc_params *params;
    jabs_stop stop; /* Stopping calculation parameters and data, set on sim_workspace_init() */
    jabs_stop stragg; /* Straggling calculation parameters and data */
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
int sim_workspace_print_spectra(const result_spectra *spectra, const char *filename);
int sim_workspace_print_bricks(const sim_workspace *ws, const char *filename);
#ifdef __cplusplus
}
#endif
#endif //JABS_SIMULATION_WORKSPACE_H
