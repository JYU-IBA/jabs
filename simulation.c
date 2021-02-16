/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#include "simulation.h"

sim_workspace *sim_workspace_init(const simulation *sim, sample *sample, jibal_gsto *gsto) {
    sim_workspace *ws = malloc(sizeof(sim_workspace));
    ws->n_reactions = sim->n_reactions;
    ws->c = calloc(sample->n_isotopes, sizeof(double));
    ws->gsto = gsto;
    ws->rk4 = 1;
    ws->stopping_type = GSTO_STO_TOT;
    ws->i_range_accel = 0;
    ws->c_x = 0.0;
    get_concs(ws, sample, ws->c_x, ws->c);
    ws->histos = calloc(ws->n_reactions, sizeof(gsl_histogram *));
    return ws;
}

void sim_workspace_free(sim_workspace *ws) {
    int i;
    for(i = 0; i < ws->n_reactions; i++) {
        gsl_histogram_free(ws->histos[i]);
    }
    free(ws->c);
    free(ws->histos);
    free(ws);
}

void sim_set_calibration(simulation *sim, double slope, double offset) {
    sim->energy_slope = slope;
    sim->energy_offset = 0.0*C_KEV;
    sim_recalculate_calibration(sim);
}

void sim_recalculate_calibration(simulation *sim) {
    sim->n_channels= ceil((1.1 * sim->ion.E - sim->energy_offset)/ sim->energy_slope);
}
