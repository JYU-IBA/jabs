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
#include "rotate.h"

simulation *sim_init() {
    simulation *sim = malloc(sizeof(simulation));
    sim->sample_theta = ALPHA; /* These defaults are for IBM geometry */
    sim->sample_phi = 0.0;
    sim->alpha = ALPHA;

    sim->theta = THETA;
    sim->beta = 180.0*C_DEG-THETA-ALPHA; /* Note: check if this is sane is defaults are changed. This should be (is) recalculated before running simulations. */
    sim->p_sr = PARTICLES_SR;
    sim->det = detector_default();
    sim->stop_step_incident = STOP_STEP_INCIDENT;
    sim->stop_step_exiting = STOP_STEP_EXITING;
    sim->fast = 0;
    sim->beam_E = ENERGY;
    sim->beam_E_broad = 0.0;

    sim->emin = E_MIN;
    sim->depthsteps_max = 0; /* Zero: automatic */
    sim->ds = 1;
    sim->ds_steps_polar = 10;
    sim->ds_steps_azi = 12;
    sim->n_ds = sim->ds_steps_polar * sim->ds_steps_azi;
    return sim;
}

void sim_free(simulation *sim) {
    free(sim);
}

void sim_calculate_geometry(simulation *sim) {
    double theta, phi; /* Temporary variables */
    rotate(0.0, 0.0, sim->sample_theta, sim->sample_phi, &theta, &phi); /* Sample in beam system. */
    sim->alpha = theta; /* SimNRA convention, alpha has no sign. */
    rotate(sim->det.theta, sim->det.phi, sim->sample_theta, sim->sample_phi, &theta, &phi); /* Detector in sample coordinate system, angles are detector in sample system. Note that for Cornell geometry phi = 90.0 deg! */
    sim->beta = C_PI - theta;
    sim->theta = sim->det.theta;
}

int sim_sanity_check(const simulation *sim) {
    if(!sim->beam_isotope) {
        fprintf(stderr, "No valid isotope given for the beam.\n");
        return -1;
    }
    if (sim->beam_E > 1000.0*C_MEV || sim->beam_E < 10*C_KEV) {
        fprintf(stderr, "Hmm...? Check your numbers. Your energy is %.5lf MeV!\n", sim->beam_E);
        return -1;
    }
    if(sim->p_sr < 0.0) {
        fprintf(stderr, "Fluence is negative (%g).\n", sim->p_sr);
        return -1;
    }
    if(detector_sanity_check(&sim->det)) {
        fprintf(stderr, "Detector failed sanity check.\n");
        return -1;
    }
    return 0;
}

#if 0
void sim_workspace_reset(sim_workspace *ws, const simulation *sim) {
    ws->i_range_accel = 0;
    ws->c_x = 0.0;
    /* TODO: finish */
}
#endif

sim_workspace *sim_workspace_init(const simulation *sim, const reaction *reactions, const sample *sample, const jibal *jibal) {
    sim_workspace *ws = malloc(sizeof(sim_workspace));
    ws->sim = *sim;
    ws->n_reactions = reaction_count(reactions);
    ws->gsto = jibal->gsto;
    ws->jibal_config = jibal->config;
    ws->isotopes = jibal->isotopes;
    ws->c_x = 0.0;

    ion_reset(&ws->ion);
    ion_set_isotope(&ws->ion, sim->beam_isotope);
    ws->ion.E = ws->sim.beam_E;
    ws->ion.S = ws->sim.beam_E_broad;

    int n_isotopes=0; /* TODO: calculate this somewhere else */
    jibal_isotope *isotope;
    for(isotope=jibal->isotopes; isotope->A != 0; isotope++) {
        n_isotopes++;
    }

    ion_nuclear_stop_fill_params(&ws->ion, jibal->isotopes, n_isotopes);

    ws->stopping_type = GSTO_STO_TOT;
    if (sim->fast) {
        ws->rk4 = 0;
        ws->nucl_stop_accurate = 0;
    } else {
        ws->rk4 = 1;
        ws->nucl_stop_accurate = 1;
    }
    sim_workspace_recalculate_calibration(ws, sim);

    if(ws->n_channels == 0) {
        free(ws);
        return NULL;
    }
    ws->reactions = calloc(ws->n_reactions, sizeof (sim_reaction));
    for(size_t i_reaction = 0; i_reaction < ws->n_reactions; i_reaction++) {
        sim_reaction *r = &ws->reactions[i_reaction];
        r->r = &reactions[i_reaction];
        ion *p = &r->p;
        ion_reset(p);
        if(sim->depthsteps_max) {
            r->n_bricks = sim->depthsteps_max;
        } else {
            if(sim->stop_step_incident == 0.0) { /* Automatic incident step size */
                r->n_bricks = (int) ceil(sim->beam_E / sqrt(sim->det.resolution) + sample->n_ranges); /* This is conservative */
            } else {
                r->n_bricks = (int) ceil(sim->beam_E / sim->stop_step_incident + sample->n_ranges); /* This is conservative */
            }
        }
        assert(r->n_bricks > 0 && r->n_bricks < 10000);
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "Number of bricks for reaction %i: %lu\n", i_reaction, r->n_bricks);
#endif
        r->histo = gsl_histogram_alloc(ws->n_channels); /* free'd by sim_workspace_free */
        gsl_histogram_set_ranges_uniform(r->histo,
                                         detector_calibrated(&sim->det, 0),
                                         detector_calibrated(&sim->det, ws->n_channels));
        r->bricks = calloc(r->n_bricks, sizeof(brick));
        if(r->r->type == REACTION_RBS) {
            ion_set_isotope(p, ws->ion.isotope);
            p->nucl_stop_isotopes = ws->ion.nucl_stop_isotopes;
            p->nucl_stop = ws->ion.nucl_stop; /* Shallow copy! Shared. */
        }
        if(r->r->type == REACTION_ERD) {
            ion_set_isotope(p, sample->isotopes[r->r->i_isotope]);
            ion_nuclear_stop_fill_params(p, jibal->isotopes, n_isotopes); /* This allocates memory */
        }
    }
    ws->c = calloc(sample->n_isotopes, sizeof(double));
    size_t i = 0;
    depth d;
    d.x = ws->c_x;
    d.i = 0;
    get_concs(sample, d, ws->c);
    return ws;
}

void sim_workspace_free(sim_workspace *ws) {
    if(!ws)
        return;
    for(size_t i = 0; i < ws->n_reactions; i++) {
        sim_reaction *r = &ws->reactions[i];
        if(r->histo) {
            gsl_histogram_free(r->histo);
            r->histo = NULL;
        }
        if(r->bricks) {
            free(r->bricks);
            r->bricks = NULL;
        }
        if(r->r->type != REACTION_RBS) {
            free(r->p.nucl_stop); /* RBS ions share nuclear stopping table. Others needs to free the memory. */
        }
    }
    free(ws->ion.nucl_stop);
    free(ws->c);
    free(ws->reactions);
    free(ws);
}

void sim_workspace_recalculate_calibration(sim_workspace *ws, const simulation *sim) {
    ws->n_channels = ceil((1.1 * sim->beam_E - sim->det.offset) / sim->det.slope);
    if(ws->n_channels > 100000) {
        ws->n_channels = 0;
    }
}

void simulation_print(FILE *f, const simulation *sim) {
    fprintf(f, "ion = %s (Z = %i, A = %i, mass %.3lf u)\n", sim->beam_isotope->name, sim->beam_isotope->Z, sim->beam_isotope->A, sim->beam_isotope->mass/C_U);
    fprintf(f, "E = %.3lf keV\n", sim->beam_E/C_KEV);
    fprintf(f, "alpha = %.3lf deg\n", sim->alpha/C_DEG);
    fprintf(f, "beta = %.3lf deg\n", sim->beta/C_DEG);
    fprintf(f, "theta = %.3lf deg\n", sim->theta/C_DEG);
    fprintf(f, "particles * sr = %e\n", sim->p_sr);
    fprintf(f, "detector calibration offset = %.3lf keV\n", sim->det.offset/C_KEV);
    fprintf(f, "detector calibration slope = %.5lf keV\n", sim->det.slope/C_KEV);
    fprintf(f, "detector resolution = %.3lf keV FWHM\n", sqrt(sim->det.resolution)*C_FWHM/C_KEV);
    fprintf(f, "step for incident ions = %.3lf keV\n", sim->stop_step_incident/C_KEV);
    fprintf(f, "step for exiting ions = %.3lf keV\n", sim->stop_step_exiting/C_KEV);
    fprintf(f, "fast level = %i\n", sim->fast);
}

void convolute_bricks(sim_workspace *ws) {
    for(size_t i = 0; i < ws->n_reactions; i++) {
        sim_reaction *r = &ws->reactions[i];
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "Reaction %i:\n", i);
#endif
        brick_int2(r->histo, r->bricks, r->n_bricks, ws->sim.det.resolution, ws->sim.p_sr);
    }
}
