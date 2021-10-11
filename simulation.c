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
#include <string.h>
#include "simulation.h"
#include "defaults.h"
#include "rotate.h"
#include "spectrum.h"
#include "message.h"

simulation *sim_init(jibal_isotope *isotopes) {
    simulation *sim = malloc(sizeof(simulation));
    sim->beam_isotope = jibal_isotope_find(isotopes, NULL, 2, 4);
    sim->sample_theta = ALPHA; /* These defaults are for IBM geometry */
    sim->sample_phi = 0.0;
    sim->fluence = FLUENCE;
    sim->beam_E = ENERGY;
    sim->emin = E_MIN;
    sim->beam_E_broad = 0.0;
    sim->channeling_offset = 1.0;
    sim->channeling_slope = 0.0;
    sim->sample = NULL; /* Needs to be set (later) */
    sim->reactions = NULL; /* Needs to be initialized after sample has been set */
    sim->n_reactions = 0;
    sim->n_det = 0;
    sim->det = NULL;
    sim->params = sim_calc_params_defaults(FALSE, FALSE);
    sim->rbs = TRUE;
    sim->erd = TRUE;
    sim_det_add(sim, detector_default(NULL));
    return sim;
}

void sim_free(simulation *sim) {
    if(!sim)
        return;
    sample_free(sim->sample);
    if(sim->det) {
        for(size_t i = 0; i < sim->n_det; i++) {
            detector_free(sim->det[i]);
        }
        free(sim->det);
    }
    for(size_t i = 0; i < sim->n_reactions; i++) {
        reaction_free(&sim->reactions[i]);
    }
    free(sim->reactions);
    free(sim);
}

sim_calc_params sim_calc_params_defaults(int ds, int fast) {
    sim_calc_params p;
    p.ds = ds;
    p.ds_steps_azi = DUAL_SCATTER_POLAR_STEPS;
    p.ds_steps_polar = DUAL_SCATTER_AZI_STEPS;
    p.n_ds = p.ds_steps_azi *  p.ds_steps_polar;
    p.stop_step_incident = STOP_STEP_INCIDENT;
    p.stop_step_exiting = STOP_STEP_EXITING;
    p.cs_n_steps = CS_CONC_STEPS;
    p.cs_stragg_half_n = CS_STRAGG_HALF_N;
    p.depthsteps_max = 0; /* automatic */
    if (fast) {
        p.rk4 = FALSE;
        p.nucl_stop_accurate = FALSE;
        p.mean_conc_and_energy = TRUE;
    } else {
        p.rk4 = TRUE;
        p.nucl_stop_accurate = TRUE;
        p.mean_conc_and_energy = FALSE;
    }
    p.cs_frac = 1.0/(1.0*(p.cs_n_steps+1));
    assert(p.cs_stragg_half_n >= 0);
    p.cs_n_stragg_steps = p.cs_stragg_half_n * 2 + 1;
    return p;
}

int sim_reactions_add_r33(simulation *sim, const jibal_isotope *jibal_isotopes, const char *filename) {
    r33_file *rfile = r33_file_read(filename);
    if(!rfile) {
        return -1;
    }
    reaction *reaction_from_file = r33_file_to_reaction(jibal_isotopes, rfile);
    if(!reaction_from_file) {
        r33_file_free(rfile);
        return -1;
    }
    fprintf(stderr, "File: %s has a reaction with %s -> %s, product %s, theta %g deg\n", filename,
            reaction_from_file->incident->name, reaction_from_file->target->name, reaction_from_file->product->name, reaction_from_file->theta/C_DEG);
    for(size_t i_reaction = 0; i_reaction < sim->n_reactions; i_reaction++) {
        reaction *r = &sim->reactions[i_reaction];
        if(reaction_is_same(r, reaction_from_file)) {
            jabs_message(MSG_INFO, stderr, "Replacing reaction %zu (%s with %s).\n", i_reaction, reaction_name(r), r->target->name);
            reaction_from_file->cs = r->cs; /* Adopt fallback cross-section from the reaction we are replacing */
            reaction_free(r);
            *r = *reaction_from_file;
            return EXIT_SUCCESS;
        }
    }
    /* Did not found a reaction to replace, just adding it. */
    sim->n_reactions++;
    sim->reactions = realloc(sim->reactions, sim->n_reactions*sizeof(reaction));
    sim->reactions[sim->n_reactions - 1] = *reaction_from_file;
    jabs_message(MSG_INFO, stderr, "Added reaction %zu.\n", sim->n_reactions - 1);
    return EXIT_SUCCESS;
}

int sim_reactions_add(simulation *sim, const sample_model *sm, reaction_type type, jibal_cross_section_type cs, double theta) { /* Note that sim->ion needs to be set! */
    if(!sim || !sim->beam_isotope || !sm) {
        return -1;
    }
    if(type == REACTION_NONE || cs ==  JIBAL_CS_NONE) {
        return 0;
    }
    struct sample *sample = sample_from_sample_model(sm);
    if(!sample) {
        return -1;
    }
#if 0
    if(type == REACTION_ERD && theta > C_PI/2.0) {
        return 0;
    }
#endif
    size_t n_reactions = sim->n_reactions + sample->n_isotopes; /* New maximum */
    sim->reactions = realloc(sim->reactions, n_reactions*sizeof(reaction));
    if(!sim->reactions) {
        sim->n_reactions = 0;
        sample_free(sample);
        return -1;
    }
    for (size_t i = 0; i < sample->n_isotopes; i++) {
        reaction *r_new = reaction_make(sim->beam_isotope, sample->isotopes[i], type, cs);
        if (!r_new) {
            jabs_message(MSG_ERROR, stderr, "Failed to make an %s reaction with isotope %zu (%s)\n", jibal_cross_section_name(cs), i, sample->isotopes[i]->name);
            continue;
        }
        sim->reactions[sim->n_reactions] = *r_new;
        free(r_new);
        sim->n_reactions++;
    };
    sample_free(sample);
    return 0;
}

void sim_reactions_free(simulation *sim) {
    if(!sim)
        return;
    for(size_t i = 0; i < sim->n_reactions; i++) {
        reaction_free(&sim->reactions[i]);
    }
    free(sim->reactions);
    sim->reactions = NULL;
    sim->n_reactions = 0;
}

int sim_sanity_check(const simulation *sim) {
    if(!sim->beam_isotope) {
        jabs_message(MSG_ERROR, stderr,  "No valid isotope given for the beam.\n");
        return -1;
    }
    if (sim->beam_E > 1000.0*C_MEV || sim->beam_E < 10*C_KEV) {
        jabs_message(MSG_ERROR, stderr,  "Hmm...? Check your numbers. Your energy is %.5lf MeV!\n", sim->beam_E);
        return -1;
    }
    if(sim->fluence < 0.0) {
        jabs_message(MSG_ERROR, stderr,  "Fluence is negative (%g).\n", sim->fluence);
        return -1;
    }
#ifdef NO_BLOCKING
    if(sim->channeling > 1.0) {
        fprintf(stderr, "Channeling correction is above 1.0 (%g)\n", sim->channeling);
        return -1;
    }
#endif
    return 0;
}

#if 0
void sim_workspace_reset(sim_workspace *ws, const simulation *sim) {
    ws->i_range_accel = 0;
    ws->c_x = 0.0;
    /* TODO: finish */
}
#endif

detector *sim_det(const simulation *sim, size_t i_det) {
    if(!sim->det)
        return NULL;
    if(i_det >= sim->n_det)
        return NULL;
    return sim->det[i_det];
}

int sim_det_add(simulation *sim, detector *det) {
    if(!sim) {
        return EXIT_FAILURE;
    }
    if(!det) {
        return EXIT_SUCCESS; /* Not an error */
    }
    sim->n_det++;
    sim->det = realloc(sim->det, sizeof(detector *) * sim->n_det);
    if(!sim->det)
        return EXIT_FAILURE;
    sim->det[sim->n_det - 1] = det;
    return EXIT_SUCCESS;
}

int sim_det_set(simulation *sim, detector *det, size_t i_det) {
    if(!sim) {
        return EXIT_FAILURE;
    }
    if(!det) {
        return EXIT_FAILURE;
    }
    if(i_det >= sim->n_det)
        return EXIT_FAILURE;
    detector_free(sim->det[i_det]);
    sim->det[i_det] = det;
    return EXIT_SUCCESS;
}

sim_workspace *sim_workspace_init(const jibal *jibal, const simulation *sim, const detector *det) {
    if(!jibal || !sim || !det) {
        jabs_message(MSG_ERROR, stderr,  "No JIBAL, sim or det. Guru thinks: %p, %p %p.\n", jibal, sim, det);
        return NULL;
    }
    if(!sim->sample) {
        jabs_message(MSG_ERROR, stderr,  "No sample has been set. Will not initialize workspace.\n");
        return NULL;
    }
    sim_workspace *ws = malloc(sizeof(sim_workspace));
    ws->sim = sim;
    ws->fluence = sim->fluence;
    ws->det = det;
    ws->sample = sim->sample;
    ws->params = sim->params;
    ws->n_reactions = sim->n_reactions;
    ws->gsto = jibal->gsto;
    ws->jibal_config = jibal->config;
    ws->isotopes = jibal->isotopes;

    if(ws->n_reactions == 0) {
        jabs_message(MSG_ERROR, stderr,  "No reactions! Will not initialize workspace if there is nothing to simulate.\n");
        free(ws);
        return NULL;
    }

    ion_reset(&ws->ion);
    ion_set_isotope(&ws->ion, sim->beam_isotope);
    ws->ion.E = ws->sim->beam_E;
    ws->ion.S = ws->sim->beam_E_broad;

    int n_isotopes=0; /* TODO: calculate this somewhere else */
    jibal_isotope *isotope;
    for(isotope=jibal->isotopes; isotope->A != 0; isotope++) {
        n_isotopes++;
    }

    ion_nuclear_stop_fill_params(&ws->ion, jibal->isotopes, n_isotopes);

    ws->stopping_type = GSTO_STO_TOT;
    sim_workspace_recalculate_n_channels(ws, sim);

    if(ws->n_channels == 0) {
        free(ws);
        jabs_message(MSG_ERROR, stderr,  "Number of channels in workspace is zero. Aborting initialization.\n");
        return NULL;
    }
    ws->histo_sum = gsl_histogram_alloc(ws->n_channels);
    if(!ws->histo_sum) {
        return NULL;
    }
    spectrum_set_calibration(ws->histo_sum, ws->det); /* Calibration can be set however already */
    gsl_histogram_reset(ws->histo_sum); /* This is not necessary, since contents should be set after simulation is over (successfully). */

    ws->reactions = calloc(ws->n_reactions, sizeof (sim_reaction));
    for(size_t i_reaction = 0; i_reaction < ws->n_reactions; i_reaction++) {
        sim_reaction *r = &ws->reactions[i_reaction];
        r->r = &sim->reactions[i_reaction];
        assert(r->r);
        assert(r->r->product);
        ion *p = &r->p;
        ion_reset(p);
        r->max_depth = 0.0;
        r->i_isotope = ws->sample->n_isotopes; /* Intentionally not valid */
        for(size_t i_isotope = 0; i_isotope < ws->sample->n_isotopes; i_isotope++) {
            if(ws->sample->isotopes[i_isotope] == r->r->target) {
#ifdef DEBUG
                fprintf(stderr, "Reaction %zu target isotope %s is isotope number %zu in sample.\n", i_reaction, r->r->target->name, i_isotope);
#endif
                r->i_isotope = i_isotope;
            }
        }
        if(ws->params.depthsteps_max) {
            r->n_bricks = ws->params.depthsteps_max;
        } else {
            if(ws->params.stop_step_incident == 0.0) { /* Automatic incident step size */
                r->n_bricks = (int) ceil(sim->beam_E / (STOP_STEP_AUTO_FUDGE_FACTOR*sqrt(ws->det->resolution)) + ws->sample->n_ranges); /* This is conservative */
            } else {
                r->n_bricks = (int) ceil(sim->beam_E / ws->params.stop_step_incident + ws->sample->n_ranges); /* This is conservative */
            }
            if(r->n_bricks > 10000) {
                jabs_message(MSG_WARNING, stderr,  "Caution: large number of bricks will be used in the simulation (%zu).\n", r->n_bricks);
            }
        }
        assert(r->n_bricks > 0 && r->n_bricks < 10000);
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "Number of bricks for reaction %i: %lu\n", i_reaction, r->n_bricks);
#endif
        r->histo = gsl_histogram_alloc(ws->n_channels); /* free'd by sim_workspace_free */
        spectrum_set_calibration(r->histo, ws->det);
        gsl_histogram_reset(r->histo);
        r->bricks = calloc(r->n_bricks, sizeof(brick));
        ion_set_isotope(p, r->r->product);
        if(r->r->type == REACTION_RBS) {
            assert(p->isotope == ws->ion.isotope);
            p->nucl_stop_isotopes = ws->ion.nucl_stop_isotopes;
            p->nucl_stop = ws->ion.nucl_stop; /* Shallow copy! Shared. */
            r->cross_section = sim_reaction_cross_section_rutherford;
        } else if(r->r->type == REACTION_ERD) {
            ion_nuclear_stop_fill_params(p, jibal->isotopes, n_isotopes); /* This allocates memory */
            r->cross_section = sim_reaction_cross_section_rutherford;
        } else if(r->r->type == REACTION_FILE) {
            ion_nuclear_stop_fill_params(p, jibal->isotopes, n_isotopes); /* This allocates memory. We could share (like with RBS) in some cases, but that's not necessarily convenient. */
            r->cross_section = sim_reaction_cross_section_tabulated;
        } else {
            p->nucl_stop = NULL;
            p->nucl_stop_isotopes = 0;
            r->cross_section = NULL;
        }
    }
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
    free(ws->reactions);
    free(ws);
}

void sim_workspace_recalculate_n_channels(sim_workspace *ws, const simulation *sim) { /* TODO: assumes calibration function is increasing */
    size_t i=0;
    while(detector_calibrated(ws->det, i) < 1.1*sim->beam_E && i <= 1000000) {i++;}
#ifdef DEBUG
    fprintf(stderr, "Simulating %zu channels\n", i);
#endif
    if(i == 1000000)
        i=0;
    ws->n_channels = i;
}

void sim_workspace_calculate_sum_spectra(sim_workspace *ws) {
    double sum;
    for(size_t i = 0; i < ws->histo_sum->n; i++) {
        sum = 0.0;
        for(size_t j = 0; j < ws->n_reactions; j++) {
            if(ws->reactions[j].histo && i < ws->reactions[j].histo->n)
                sum += ws->reactions[j].histo->bin[i];
        }
        ws->histo_sum->bin[i] = sum;
    }
}

void simulation_print(FILE *f, const simulation *sim) {
    if(!sim) {
        return;
    }
    double theta, phi; /* Temporary variables */
    double alpha; /* Incident angle, SimNRA convention (no signs). */
    rotate(0.0, 0.0, sim->sample_theta, sim->sample_phi, &theta, &phi); /* Sample in beam system. */
    alpha = theta;
    if(sim->beam_isotope) {
        jabs_message(MSG_INFO, f,  "ion = %s (Z = %i, A = %i, mass %.3lf u)\n", sim->beam_isotope->name, sim->beam_isotope->Z, sim->beam_isotope->A, sim->beam_isotope->mass / C_U);
    } else {
        jabs_message(MSG_INFO, f, "ion = None\n");
    }
    jabs_message(MSG_INFO, f, "E = %.3lf keV\n", sim->beam_E/C_KEV);
    jabs_message(MSG_INFO, f, "E_broad = %.3lf keV FWHM\n", sqrt(sim->beam_E_broad)*C_FWHM/C_KEV);
    jabs_message(MSG_INFO, f, "E_min = %.3lf keV\n", sim->emin/C_KEV);
    jabs_message(MSG_INFO, f, "alpha = %.3lf deg\n", alpha/C_DEG);
    jabs_message(MSG_INFO, f, "n_detectors = %zu\n", sim->n_det);
    jabs_message(MSG_INFO, f, "n_reactions = %zu\n", sim->n_reactions);
    jabs_message(MSG_INFO, f, "fluence = %e (%.5lf p-uC)\n", sim->fluence, sim->fluence*C_E*1.0e6);
    jabs_message(MSG_INFO, f, "step for incident ions = %.3lf keV\n", sim->params.stop_step_incident/C_KEV);
    jabs_message(MSG_INFO, f, "step for exiting ions = %.3lf keV\n", sim->params.stop_step_exiting/C_KEV);
    jabs_message(MSG_INFO, f, "stopping RK4 = %s\n", sim->params.rk4?"true":"false");
    jabs_message(MSG_INFO, f, "depth steps max = %zu\n", sim->params.depthsteps_max);
    jabs_message(MSG_INFO, f, "cross section weighted by straggling = %s\n", sim->params.mean_conc_and_energy?"false":"true");
    jabs_message(MSG_INFO, f, "accurate nuclear stopping = %s\n", sim->params.nucl_stop_accurate?"true":"false");
    if(sim->channeling_offset != 1.0 || sim->channeling_slope != 0.0) {
        jabs_message(MSG_INFO, f, "substrate channeling yield correction offset = %.5lf\n", sim->channeling_offset);
        jabs_message(MSG_INFO, f, "substrate channeling yield correction slope = %g / keV (%e)\n", sim->channeling_slope/(1.0/C_KEV), sim->channeling_slope);
    }
}

void convolute_bricks(sim_workspace *ws) {
    for(size_t i = 0; i < ws->n_reactions; i++) {
        sim_reaction *r = &ws->reactions[i];
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "Reaction %i:\n", i);
#endif
        brick_int2(r->histo, r->bricks, r->last_brick, ws->det->resolution, ws->fluence * ws->det->solid);
    }
}

void sim_reaction_recalculate_internal_variables(sim_reaction *sim_r) { /* Calculate variables for Rutherford (and Andersen) cross sections. This is done for all reactions. */
    const jibal_isotope *incident = sim_r->r->incident;
    const jibal_isotope *target = sim_r->r->target;
    const jibal_isotope *product = sim_r->r->product;
    sim_r->E_cm_ratio = target->mass / (incident->mass + target->mass);
    sim_r->mass_ratio = incident->mass / target->mass;
    if(product == incident) { /* RBS */
        if(incident->mass >= target->mass && sim_r->theta > asin(target->mass / incident->mass)) {
            sim_r->K = 0.0;
            sim_r->cs_constant = 0.0;
            sim_r->stop = TRUE;
            return;
        }
        sim_r->K = jibal_kin_rbs(sim_r->r->incident->mass, sim_r->r->target->mass, sim_r->theta, '+');
        sim_r->theta_cm = sim_r->theta + asin(sim_r->mass_ratio * sin(sim_r->theta));
        sim_r->cs_constant = (pow2(sin(sim_r->theta_cm))) / (pow2(sin(sim_r->theta)) * cos(sim_r->theta_cm - sim_r->theta)) *
                             pow2((incident->Z * C_E * target->Z * C_E) / (4.0 * C_PI * C_EPSILON0)) *
                             pow4(1.0 / sin(sim_r->theta_cm / 2.0)) * (1.0 / 16.0);
    } else if(product == target) { /* ERD */
        if(sim_r->theta > C_PI/2.0) {
#ifdef DEBUG
            fprintf(stderr, "ERD with %s is not possible (theta %g deg > 90.0 deg)\n", target->name, sim_r->theta);
#endif
            sim_r->K = 0.0;
            sim_r->cs_constant = 0.0;
            sim_r->theta_cm = 0.0;
            sim_r->stop = TRUE;
            return;
        }
        sim_r->K = jibal_kin_erd(sim_r->r->incident->mass, sim_r->r->target->mass, sim_r->theta);
        sim_r->theta_cm = C_PI - 2.0 * sim_r->theta;
        sim_r->cs_constant = pow2(incident->Z * C_E * target->Z * C_E / (8 * C_PI * C_EPSILON0)) * pow2(1.0 + incident->mass / target->mass) * pow(cos(sim_r->theta), -3.0) * pow2(sim_r->E_cm_ratio);
    } else {
        sim_r->K = 0.0;
        sim_r->cs_constant = 0.0;
        sim_r->theta_cm = 0.0;
        sim_r->stop = TRUE;
    }

    sim_r->r_VE_factor = 48.73 * C_EV * incident->Z * target->Z * sqrt(pow(incident->Z, 2.0 / 3.0) + pow(target->Z, 2.0 / 3.0)); /* Factors for Andersen correction */
    sim_r->r_VE_factor2 = pow2(0.5 / sin(sim_r->theta_cm / 2.0));
#ifdef DEBUG
    fprintf(stderr, "Reaction recalculated, theta = %g deg, theta_cm = %g deg, K = %g\n", sim_r->theta/C_DEG, sim_r->theta_cm/C_DEG, sim_r->K);
#endif
}

void sim_reaction_reset_bricks(sim_reaction *sim_r) {
    memset(sim_r->bricks, 0, sizeof(brick) * sim_r->n_bricks);
}

double sim_reaction_andersen(const sim_reaction *sim_r, double E_cm) {
    const double r_VE = sim_r->r_VE_factor / E_cm;
    return pow2(1 + 0.5 * r_VE) / pow2(1 + r_VE + sim_r->r_VE_factor2 * pow2(r_VE));
}

double sim_reaction_cross_section_rutherford(const sim_reaction *sim_r, double E) {
#ifdef CROSS_SECTIONS_FROM_JIBAL
    return jibal_cross_section_erd(sim_r->r->incident, sim_r->r->target, sim_r->theta, E, sim_r->r->cs);
#else
    const reaction *r = sim_r->r;
    const double E_cm = sim_r->E_cm_ratio * E;
    double sigma_r = sim_r->cs_constant / pow2(E_cm) ;
    switch(r->cs) {
        case JIBAL_CS_RUTHERFORD:
            return sigma_r;
        case JIBAL_CS_ANDERSEN:
            return sigma_r * sim_reaction_andersen(sim_r, E_cm);
        default:
            return 0.0;
    }
#endif
}

double sim_reaction_cross_section_tabulated(const sim_reaction *sim_r, double E) {
    size_t lo, mi, hi;
    const reaction *r = sim_r->r;
    const struct reaction_point *t = r->cs_table;
    hi = r->n_cs_table - 1;
    lo = 0;
    if(E < t[lo].E || E > t[hi].E) {
        return sim_reaction_cross_section_rutherford(sim_r, E); /* Fall back quietly to analytical formulae outside tabulated values */
    }
    if(sim_r->theta != sim_r->r->theta) {
#ifdef DEBUG_VERBOSE /* This happens with DS */
        fprintf(stderr, "Reaction theta %g deg different from one in the file: %g deg\n", sim_r->theta/C_DEG, sim_r->r->theta/C_DEG);
#endif
        return sim_reaction_cross_section_rutherford(sim_r, E); /* Fall back quietly if reaction theta has been changed from original scattering angle (in the file) */
    }
    while (hi - lo > 1) {
        mi = (hi + lo) / 2;
        if (E >= t[mi].E) {
            lo = mi;
        } else {
            hi = mi;
        }
    }
    return t[lo].sigma+((t[lo+1].sigma-t[lo].sigma)/(t[lo+1].E-t[lo].E))*(E-t[lo].E);
}

double sim_calculate_exit_angle(const simulation *sim, const detector *det) {
    double theta, phi;
    rotate(det->theta, det->phi, sim->sample_theta, sim->sample_phi, &theta, &phi); /* Detector in sample coordinate system, angles are detector in sample system. Note that for Cornell geometry phi = 90.0 deg! */
    return C_PI - theta;
}
