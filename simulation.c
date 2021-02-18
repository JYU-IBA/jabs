/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#include <assert.h>
#include "simulation.h"
#include "defaults.h"

simulation *sim_init() {
    simulation *sim = malloc(sizeof(simulation));
    sim->alpha = ALPHA;
    sim->beta = BETA;
    sim->theta = THETA;
    sim->p_sr = PARTICLES_SR;
    sim->energy_resolution = (DETECTOR_RESOLUTION*DETECTOR_RESOLUTION);
    sim->stop_step_incident = STOP_STEP_INCIDENT;
    sim->stop_step_exiting = STOP_STEP_EXITING;
    sim->fast = 0;
    sim->ion.E = ENERGY;
    sim->ion.S = 0.0;
    sim->energy_slope = ENERGY_SLOPE;
    sim->energy_offset = 0.0*C_KEV;
    sim->emin = E_MIN;
    ion_set_angle(&sim->ion, sim->alpha);
    return sim;
}

void sim_free(simulation *sim) {
    free(sim);
}

int sim_sanity_check(const simulation *sim) {
    if(!sim->ion.isotope) {
        fprintf(stderr, "No valid isotope given for the beam.\n");
        return -1;
    }
    if (sim->ion.E > 1000.0*C_MEV || sim->ion.E < 10*C_KEV) {
        fprintf(stderr, "Hmm...? Check your numbers. Your energy is %.5lf MeV!\n", sim->ion.E);
        return -1;
    }
    return 0;
}

sim_workspace *sim_workspace_init(const simulation *sim, const sample *sample, jibal_gsto *gsto) {
    sim_workspace *ws = malloc(sizeof(sim_workspace));
    ws->n_reactions = sim->n_reactions;
    ws->gsto = gsto;
    ws->i_range_accel = 0;
    ws->c_x = 0.0;

    ws->p_sr_cos_alpha = sim->p_sr / cos(sim->alpha);
    if (sim->fast) {
        ws->stopping_type = GSTO_STO_ELE;
        ws->rk4 = 0;
    } else {
        ws->stopping_type = GSTO_STO_TOT;
        ws->rk4 = 1;
    }
    sim_workspace_recalculate_calibration(ws, sim);
    if(ws->n_channels == 0) {
        free(ws);
        return NULL;
    }
    ws->histos = calloc(ws->n_reactions, sizeof(gsl_histogram *));
    int i_reaction;
    for(i_reaction = 0; i_reaction < ws->n_reactions; i_reaction++) {
        ws->histos[i_reaction] = gsl_histogram_alloc(ws->n_channels); /* free'd by sim_workspace_free */
        gsl_histogram_set_ranges_uniform(ws->histos[i_reaction],
                                         sim->energy_offset,
                                         sim->energy_offset + sim->energy_slope * ws->n_channels);
    }
    ws->c = calloc(sample->n_isotopes, sizeof(double));
    get_concs(ws, sample, ws->c_x, ws->c);
    return ws;
}

void sim_workspace_free(sim_workspace *ws) {
    if(!ws)
        return;
    int i;
    for(i = 0; i < ws->n_reactions; i++) {
        if(ws->histos[i]) {
            gsl_histogram_free(ws->histos[i]);
            ws->histos[i] = NULL;
        }
    }
    free(ws->c);
    free(ws->histos);
    free(ws);
}

void sim_set_calibration(simulation *sim, double slope, double offset) {
    sim->energy_slope = slope;
    sim->energy_offset = 0.0*C_KEV;
}

void sim_workspace_recalculate_calibration(sim_workspace *ws, const simulation *sim) {
    ws->n_channels = ceil((1.1 * sim->ion.E - sim->energy_offset) / sim->energy_slope);
    if(ws->n_channels < 0 || ws->n_channels > 100000) {
        ws->n_channels = 0;
    }
}

void simulation_print(FILE *f, const simulation *sim) {
    fprintf(stderr, "ion = %s (Z = %i, A = %i, mass %.3lf u)\n", sim->ion.isotope->name, sim->ion.isotope->Z, sim->ion.isotope->A, sim->ion.isotope->mass/C_U);
    fprintf(stderr, "E = %.3lf keV\n", sim->ion.E/C_KEV);
    fprintf(stderr, "alpha = %.3lf deg\n", sim->alpha/C_DEG);
    fprintf(stderr, "beta = %.3lf deg\n", sim->beta/C_DEG);
    fprintf(stderr, "theta = %.3lf deg\n", sim->theta/C_DEG);
    fprintf(stderr, "particles * sr = %e\n", sim->p_sr);
    fprintf(stderr, "calibration offset = %.3lf keV\n", sim->energy_offset/C_KEV);
    fprintf(stderr, "calibration slope = %.5lf keV\n", sim->energy_slope/C_KEV);
    fprintf(stderr, "detector resolution = %.3lf keV FWHM\n", sqrt(sim->energy_resolution)*C_FWHM/C_KEV);
    fprintf(stderr, "step for incident ions = %.3lf keV\n", sim->stop_step_incident/C_KEV);
    fprintf(stderr, "step for exiting ions = %.3lf keV\n", sim->stop_step_exiting/C_KEV);
    fprintf(stderr, "fast level = %i\n", sim->fast);
}
