/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

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

simulation *sim_init(jibal *jibal) {
    simulation *sim = malloc(sizeof(simulation));
    sim->beam_isotope = jibal_isotope_find(jibal->isotopes, NULL, 2, 4);
    sim->beam_aperture = NULL;
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
    sim->params = sim_calc_params_defaults(NULL);
    sim->rbs = TRUE;
    sim->erd = TRUE;
    sim->cs_rbs = jabs_reaction_cs_from_jibal_cs(jibal->config->cs_rbs);
    sim->cs_erd =  jabs_reaction_cs_from_jibal_cs(jibal->config->cs_erd);
    sim_det_add(sim, detector_default(NULL));
    return sim;
}

void sim_free(simulation *sim) {
    if(!sim)
        return;
    sample_free(sim->sample);
    aperture_free(sim->beam_aperture);
    if(sim->det) {
        for(size_t i = 0; i < sim->n_det; i++) {
            detector_free(sim->det[i]);
        }
        free(sim->det);
    }
    for(size_t i = 0; i < sim->n_reactions; i++) {
        reaction_free(sim->reactions[i]);
    }
    free(sim->reactions);
    sim_calc_params_free(sim->params);
    free(sim);
}

sim_calc_params *sim_calc_params_defaults(sim_calc_params *p) {
    if(!p) {
        p = malloc(sizeof(sim_calc_params));
        p->cs_stragg_pd = NULL; /* Will be allocated by sim_calc_params_update() */
    }
    p->stop_step_incident = STOP_STEP_INCIDENT;
    p->stop_step_exiting = STOP_STEP_EXITING;
    p->stop_step_fudge_factor = STOP_STEP_FUDGE_FACTOR;
    p->stop_step_min = STOP_STEP_MIN;
    p->stop_step_add = STOP_STEP_ADD;
    p->depthsteps_max = 0; /* automatic */
    p->geostragg = FALSE;
    p->beta_manual  = FALSE;
    p->ds = FALSE;
    p->ds_steps_azi = 0;
    p->ds_steps_polar = 0;
    p->rk4 = TRUE;
    p->nucl_stop_accurate = TRUE;
    p->mean_conc_and_energy = FALSE;
    p->cs_n_stragg_steps = CS_STRAGG_STEPS;
    p->cs_n_steps = CS_CONC_STEPS;
    p->rough_layer_multiplier = 1.0;
    p->sigmas_cutoff = SIGMAS_CUTOFF;
    p->gaussian_accurate = FALSE;
    p->int_cs_max_intervals = CS_CONC_MAX_INTEGRATION_INTERVALS;
    p->int_cs_accuracy = CS_CONC_INTEGRATION_ACCURACY;
    p->int_cs_stragg_max_intervals = CS_STRAGG_MAX_INTEGRATION_INTERVALS;
    p->int_cs_stragg_accuracy = CS_STRAGG_INTEGRATION_ACCURACY;
#ifdef DEBUG
    fprintf(stderr, "New calc params created.\n");
#endif
    return p;
}

sim_calc_params *sim_calc_params_defaults_fast(sim_calc_params *p) {
    sim_calc_params_defaults(p);
    p->rk4 = FALSE;
    p->nucl_stop_accurate = FALSE;
    p->mean_conc_and_energy = TRUE;
    p->cs_n_stragg_steps = 0; /* Not used if mean_conc_and_energy == TRUE */
    p->cs_n_steps = 0; /* Not used if mean_conc_and_energy == TRUE */
    p->stop_step_fudge_factor *= 1.4;
    p->stop_step_add *= 2.0;
    p->geostragg = FALSE;
    p->rough_layer_multiplier = 0.5;
    p->sigmas_cutoff = SIGMAS_FAST_CUTOFF;
    return p;
}

sim_calc_params *sim_calc_params_defaults_accurate(sim_calc_params *p) {
    sim_calc_params_defaults(p);
    p->cs_n_steps = 0; /* Automatic (adaptive) */
    p->cs_n_stragg_steps = 0; /* Automatic (adaptive) */
    p->stop_step_fudge_factor *= 0.5;
    p->sigmas_cutoff += 1.0;
    p->gaussian_accurate = TRUE;
    return p;
}

sim_calc_params *sim_calc_params_defaults_brisk(sim_calc_params *p) {
    sim_calc_params_defaults(p);
    p->cs_n_stragg_steps -= 2;
    p->stop_step_fudge_factor *= 1.25;
    p->stop_step_add *= 2.0;
    p->sigmas_cutoff -= 1.0;
    return p;
}

sim_calc_params *sim_calc_params_defaults_improved(sim_calc_params *p) {
    sim_calc_params_defaults(p);
    p->cs_n_steps += 2;
    p->cs_n_stragg_steps += 4;
    p->stop_step_fudge_factor *= 0.75;
    p->sigmas_cutoff += 0.5;
    p->gaussian_accurate = TRUE;
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
    p->cs_frac = 1.0/(1.0*(p->cs_n_steps+1));
    prob_dist_free(p->cs_stragg_pd);
    p->cs_stragg_pd = prob_dist_gaussian(p->cs_n_stragg_steps);

    if(p->ds && p->ds_steps_azi == 0 && p->ds_steps_polar == 0) { /* DS defaults are applied if nothing else is specified */
        p->ds_steps_azi = DUAL_SCATTER_AZI_STEPS;
        p->ds_steps_polar = DUAL_SCATTER_POLAR_STEPS;
    }
    p->n_ds = p->ds_steps_azi *  p->ds_steps_polar;
}

void sim_calc_params_ds(sim_calc_params *p, int ds) {
    if(ds) {
        p->ds = TRUE;
    }
}

void sim_calc_params_print(const sim_calc_params *params) {
    if(!params)
        return;
    jabs_message(MSG_INFO, stderr, "step for incident ions = %.3lf keV (0 = auto)\n", params->stop_step_incident/C_KEV);
    jabs_message(MSG_INFO, stderr, "step for exiting ions = %.3lf keV (0 = auto)\n", params->stop_step_exiting/C_KEV);
    jabs_message(MSG_INFO, stderr, "stopping step fudge factor = %g\n", params->stop_step_fudge_factor);
    jabs_message(MSG_INFO, stderr, "stopping step minimum = %.3lf keV (0 = auto)\n", params->stop_step_min / C_KEV);
    jabs_message(MSG_INFO, stderr, "stopping RK4 = %s\n", params->rk4?"true":"false");
    jabs_message(MSG_INFO, stderr, "accurate nuclear stopping = %s\n", params->nucl_stop_accurate?"true":"false");
    jabs_message(MSG_INFO, stderr, "depth steps max = %zu\n", params->depthsteps_max);
    jabs_message(MSG_INFO, stderr, "cross section of brick determined using mean concentration and energy = %s\n", params->mean_conc_and_energy?"true":"false");
    jabs_message(MSG_INFO, stderr, "geometric broadening = %s\n", params->geostragg?"true":"false");
    if(!params->mean_conc_and_energy) {
        if(params->cs_n_steps == 0) {
            jabs_message(MSG_INFO, stderr, "concentration * cross section steps integration accuracy = %g\n", params->int_cs_accuracy);
        } else {
            jabs_message(MSG_INFO, stderr, "concentration * cross section steps per brick = %zu\n", params->cs_n_steps);
        }
        if(params->cs_n_stragg_steps == 0) {
            jabs_message(MSG_INFO, stderr, "straggling weighting integration accuracy = %g\n", params->int_cs_stragg_accuracy);
        } else {
            jabs_message(MSG_INFO, stderr, "straggling substeps = %zu\n", params->cs_n_stragg_steps);
        }
    }
}

jabs_reaction_cs sim_cs(const simulation *sim, const reaction_type type) {
    if(type == REACTION_RBS ||type == REACTION_RBS_ALT)
        return sim->cs_rbs;
    if(type == REACTION_ERD)
        return sim->cs_erd;
    return JABS_CS_NONE;
}

int sim_reactions_add_reaction(simulation *sim, reaction *r) {
    if(!sim || !r)
        return EXIT_FAILURE;
    sim->n_reactions++;
    sim->reactions = realloc(sim->reactions, sim->n_reactions*sizeof(reaction *));
    sim->reactions[sim->n_reactions - 1] = r;
    jabs_message(MSG_INFO, stderr, "Added reaction %zu (%s), %s(%s,%s)%s, E = [%g keV, %g keV], Q = %g keV\n",
                 sim->n_reactions, reaction_name(r), r->target->name, r->incident->name, r->product->name, r->product_nucleus->name, r->E_min / C_KEV, r->E_max / C_KEV, r->Q / C_KEV);
    return EXIT_SUCCESS;
}

int sim_reactions_remove_reaction(simulation *sim, size_t i) {
    if(i >= sim->n_reactions) {
        jabs_message(MSG_ERROR, stderr, "There are %zu reactions, can't remove reaction %zu\n", sim->n_reactions, i + 1);
        return EXIT_FAILURE;
    }
    reaction_free(sim->reactions[i]);
    sim->reactions[i] = NULL;
    sim_sort_reactions(sim);
    sim->n_reactions--;
    return EXIT_SUCCESS;
}

int sim_reactions_add_r33(simulation *sim, const jibal_isotope *jibal_isotopes, const char *filename) {
    r33_file *rfile = r33_file_read(filename);
    if(!rfile) {
        jabs_message(MSG_ERROR, stderr, "Could not load R33 from file \"%s\".\n", filename);
        return EXIT_FAILURE;
    }
    reaction *r = r33_file_to_reaction(jibal_isotopes, rfile);
    r33_file_free(rfile);
    if(!r) {
        jabs_message(MSG_ERROR, stderr, "Could not convert R33 file to a reaction!\n");
        return EXIT_FAILURE;
    }
    jabs_message(MSG_INFO, stderr, "File: %s has a reaction with %s -> %s, product %s, theta %g deg\n",
            filename, r->incident->name, r->target->name, r->product->name, r->theta/C_DEG);
    sim_reactions_add_reaction(sim, r);
    return EXIT_SUCCESS;
}

int sim_reactions_add_auto(simulation *sim, const sample_model *sm, reaction_type type, jabs_reaction_cs cs) { /* Note that sim->ion needs to be set! */
    if(!sim || !sim->beam_isotope || !sm) {
        return -1;
    }
    if(type != REACTION_RBS && type != REACTION_RBS_ALT && type != REACTION_ERD) {
        return 0;
    }
    if(cs == JIBAL_CS_NONE) {
        return 0;
    }
    struct sample *sample = sample_from_sample_model(sm);
    if(!sample) {
        return -1;
    }
    for (size_t i = 0; i < sample->n_isotopes; i++) {
        const jibal_isotope *isotope = sample->isotopes[i];
        if(type == REACTION_RBS_ALT) {
            if(sim->beam_isotope->mass <= isotope->mass) { /* alternative solution for RBS exists only if m1 > m2 */
                continue;
            }
        }
        reaction *r_new = reaction_make(sim->beam_isotope, isotope, type, cs);
        if (!r_new) {
            jabs_message(MSG_ERROR, stderr, "Failed to make an %s reaction with %s cross sections for isotope %zu (%s)\n", reaction_type_to_string(type), jabs_reaction_cs_to_string(cs), i, isotope->name);
            continue;
        }
        sim_reactions_add_reaction(sim, r_new);
    };
    sample_free(sample);
    return 0;
}

void sim_reactions_free(simulation *sim) {
    if(!sim)
        return;
    for(size_t i = 0; i < sim->n_reactions; i++) {
        reaction_free(sim->reactions[i]);
    }
    free(sim->reactions);
    sim->reactions = NULL;
    sim->n_reactions = 0;
}

int sim_sanity_check(const simulation *sim) { /* This does not guarantee sanity, but it guarantees insanity... Checks should be limited to obvious things. Note that sample and reactions might not be set before this is called. */
    if(!sim) {
        jabs_message(MSG_ERROR, stderr,  "No simulation.\n");
        return -1;
    }

    if(!sim->beam_isotope) {
        jabs_message(MSG_ERROR, stderr,  "No valid isotope given for the beam.\n");
        return -1;
    }
    if(sim->beam_isotope->Z == 0) {
        jabs_message(MSG_ERROR, stderr,  "Neutron beams are not supported.\n");
        return -1;
    }
    if (sim->beam_E > 1000.0*C_MEV || sim->beam_E < 10*C_KEV) {
        jabs_message(MSG_ERROR, stderr,  "Hmm...? Check your numbers. Your energy is %g J (%g MeV)!\n", sim->beam_E, sim->beam_E/C_MEV);
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

detector *sim_det_from_string(const simulation *sim, const char *s) {
    size_t i_det;
    if(*s >= '0' && *s <= '9') {
        i_det = strtoull(s, NULL, 10);
        if(i_det == 0)
            return NULL;
        i_det--;
        return sim_det(sim, i_det);
    }
    /* TODO: other detector naming besides numbering */
    return NULL;
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

void sim_workspace_init_reactions(sim_workspace *ws) {
    const simulation *sim = ws->sim;
    ws->reactions = calloc(sim->n_reactions, sizeof (sim_reaction *));
    for(size_t i_reaction = 0; i_reaction < sim->n_reactions; i_reaction++) {
        ws->reactions[i_reaction] = sim_reaction_init(ws, sim->reactions[i_reaction]);
        ws->n_reactions++;
    }
}

void sim_workspace_calculate_number_of_bricks(sim_workspace *ws) {
    size_t n_bricks = 0;
    const detector *det = ws->det;
    const simulation *sim = ws->sim;
    if(ws->params->depthsteps_max) {
        n_bricks = ws->params->depthsteps_max;
    } else {
        if(det->type == DETECTOR_ENERGY) {
            if(ws->params->stop_step_incident == 0.0) { /* Automatic incident step size */
                n_bricks = (int) ceil(sim->beam_E / (ws->params->stop_step_fudge_factor*sqrt(detector_resolution(ws->det, sim->beam_isotope, sim->beam_E))) + ws->sample->n_ranges);
            } else {
                n_bricks = (int) ceil(sim->beam_E / ws->params->stop_step_incident + ws->sample->n_ranges); /* This is conservative */
                fprintf(stderr, "n_bricks = %zu\n", n_bricks);
            }
        } else { /* TODO: maybe something more clever is needed here */
            n_bricks = BRICKS_DEFAULT;
        }
        if(n_bricks > BRICKS_MAX) {
            n_bricks = BRICKS_MAX;
            jabs_message(MSG_WARNING, stderr,  "Number of bricks limited to %zu.\n", n_bricks);
        }
    }
    ws->n_bricks = n_bricks;
#ifdef DEBUG
    fprintf(stderr, "Number of bricks: %zu\n", n_bricks);
#endif
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
    ws->emin = sim->emin;
    ws->fluence = sim->fluence;
    ws->det = det;
    ws->sample = sim->sample;
    ws->params = sim_calc_params_defaults(NULL);
    sim_calc_params_copy(sim->params,  ws->params);
    sim_calc_params_update(ws->params);
    if(ws->params->stop_step_min == 0.0) { /* Automatic, calculate here. */
        if(det->type == DETECTOR_ENERGY) {
            ws->params->stop_step_min = sqrt(det->calibration->resolution_variance)/2.0;
        } else {
            ws->params->stop_step_min = STOP_STEP_MIN_FALLBACK;
        }
    }
#ifdef DEBUG
    fprintf(stderr, "Minimum stop step %g keV.\n", ws->params->stop_step_min/C_KEV);
#endif
    ws->gsto = jibal->gsto;
    ws->isotopes = jibal->isotopes;
    ws->n_reactions = 0; /* Will be incremented later */

    if(sim->n_reactions == 0) {
        jabs_message(MSG_ERROR, stderr,  "No reactions! Will not initialize workspace if there is nothing to simulate.\n");
        free(ws);
        return NULL;
    }

    if(sim->params->beta_manual && sim->params->ds) {
        jabs_message(MSG_WARNING, stderr,  "Manual exit angle is enabled in addition to dual scattering. This is an unsupported combination. Manual exit angle calculation will be disabled.\n");
        ws->params->beta_manual  = FALSE;
    }

    if(sim->params->beta_manual && sim->params->geostragg) {
        jabs_message(MSG_WARNING, stderr,  "Manual exit angle is enabled in addition to geometric scattering. This is an unsupported combination. Geometric straggling calculation will be disabled.\n");
        ws->params->geostragg = FALSE;
    }
    ws->ion = sim->ion; /* Shallow copy, but that is ok */

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
    spectrum_set_calibration(ws->histo_sum, ws->det, JIBAL_ANY_Z); /* Calibration (assuming default calibration) can be set now. */
    gsl_histogram_reset(ws->histo_sum); /* This is not necessary, since contents should be set after simulation is over (successfully). */

    sim_workspace_calculate_number_of_bricks(ws);
    sim_workspace_init_reactions(ws);

    if(ws->params->cs_n_steps == 0) { /* Actually integrate, allocate workspace for this */
#ifdef DEBUG
        fprintf(stderr, "cs_n_steps = 0, allocating integration workspace w_int_cs with %zu max intervals.\n", ws->params->int_cs_max_intervals);
#endif
        ws->w_int_cs = gsl_integration_workspace_alloc(ws->params->int_cs_max_intervals);
    } else {
        ws->w_int_cs = NULL;
    }
    if(ws->params->cs_n_stragg_steps == 0) {
#ifdef DEBUG
        fprintf(stderr, "cs_n_stragg_steps = 0, allocating integration workspace w_int_cs_stragg with %zu max intervals.\n", ws->params->int_cs_stragg_max_intervals);
#endif
        ws->w_int_cs_stragg = gsl_integration_workspace_alloc(ws->params->int_cs_stragg_max_intervals);
    } else {
        ws->w_int_cs_stragg = NULL;
    }
    return ws;
}

void sim_workspace_free(sim_workspace *ws) {
    if(!ws)
        return;
    sim_calc_params_free(ws->params);
    for(size_t i = 0; i < ws->n_reactions; i++) {
        sim_reaction_free(ws->reactions[i]);
    }
    free(ws->reactions);
    gsl_integration_workspace_free(ws->w_int_cs);
    gsl_integration_workspace_free(ws->w_int_cs_stragg);
    gsl_histogram_free(ws->histo_sum);
    free(ws);
}

void sim_workspace_recalculate_n_channels(sim_workspace *ws, const simulation *sim) { /* TODO: assumes calibration function is increasing */
    size_t n_max = CHANNELS_ABSOLUTE_MIN; /* Always simulate at least CHANNELS_ABSOLUTE_MIN channels */
    for(size_t i_reaction = 0; i_reaction < sim->n_reactions; i_reaction++) {
        const reaction *r = sim->reactions[i_reaction];
        if(!reaction_is_possible(r, ws->det->theta)) {
#ifdef DEBUG
            fprintf(stderr, "Reaction %zu (target %s, type %s) is not possible when theta = %g deg. Skipping. \n", i_reaction + 1, r->target->name,
                    reaction_type_to_string(r->type), ws->det->theta / C_DEG);
#endif
            continue;
        }
        double E = reaction_product_energy(r, ws->det->theta, sim->beam_E);
        double E_safer = E + 3.0*sqrt(detector_resolution(ws->det, r->product, E)); /* Add 3x resolution sigma to max energy */
        E_safer *= 1.1; /* and 10% for good measure! */
        while(detector_calibrated(ws->det, JIBAL_ANY_Z, n_max) < E_safer && n_max <= CHANNELS_ABSOLUTE_MAX) {n_max++;} /* Increase number of channels until we hit this energy. TODO: this requires changes for ToF spectra. */
#ifdef DEBUG
        fprintf(stderr, "Reaction %zu (target %s, type %s), max E = %g keV -> %g keV after resolution and safety factor have been added, current n_max %zu (channels)\n", i_reaction + 1, r->target->name,
                reaction_type_to_string(r->type), E/C_KEV, E_safer/C_KEV, n_max);
#endif
    }
    if(n_max == CHANNELS_ABSOLUTE_MAX)
        n_max=0;
    ws->n_channels = n_max;
}

void sim_workspace_calculate_sum_spectra(sim_workspace *ws) {
    double sum;
    for(size_t i = 0; i < ws->histo_sum->n; i++) {
        sum = 0.0;
        for(size_t j = 0; j < ws->n_reactions; j++) {
            if(ws->reactions[j]->histo && i < ws->reactions[j]->histo->n)
                sum += ws->reactions[j]->histo->bin[i];
        }
        ws->histo_sum->bin[i] = sum;
    }
}

void sim_print(const simulation *sim) {
    if(!sim) {
        return;
    }
    if(sim->beam_isotope) {
        jabs_message(MSG_INFO, stderr,  "ion = %s (Z = %i, A = %i, mass %.3lf u)\n", sim->beam_isotope->name, sim->beam_isotope->Z, sim->beam_isotope->A, sim->beam_isotope->mass / C_U);
    } else {
        jabs_message(MSG_INFO, stderr, "ion = None\n");
    }
    jabs_message(MSG_INFO, stderr, "E = %.3lf keV\n", sim->beam_E/C_KEV);
    jabs_message(MSG_INFO, stderr, "E_broad = %.3lf keV FWHM\n", sqrt(sim->beam_E_broad)*C_FWHM/C_KEV);
    jabs_message(MSG_INFO, stderr, "E_min = %.3lf keV\n", sim->emin/C_KEV);
    jabs_message(MSG_INFO, stderr, "alpha = %.3lf deg\n", sim_alpha_angle(sim)/C_DEG);
    jabs_message(MSG_INFO, stderr, "sample tilt (horizontal) = %.3lf deg\n", angle_tilt(sim->sample_theta, sim->sample_phi, 'x')/C_DEG);
    jabs_message(MSG_INFO, stderr, "sample tilt (vertical) = %.3lf deg\n", angle_tilt(sim->sample_theta, sim->sample_phi, 'y')/C_DEG);
    rot_vect v = rot_vect_from_angles(C_PI - sim->sample_theta, sim->sample_phi); /* By default our sample faces the beam and tilt angles are based on that choice. Pi is there for a reason. */
    jabs_message(MSG_INFO, stderr, "surf normal unit vector (beam in z direction) = (%.3lf, %.3lf, %.3lf)\n", v.x, v.y, v.z);
    char *aperture_str = aperture_to_string(sim->beam_aperture);
    jabs_message(MSG_INFO, stderr, "aperture = %s\n", aperture_str);
    free(aperture_str);
    jabs_message(MSG_INFO, stderr, "n_detectors = %zu\n", sim->n_det);
    for(size_t i = 0; i < sim->n_det; i++) {
        detector *det = sim->det[i];
        jabs_message(MSG_INFO, stderr, "DETECTOR %zu (run 'show detector %zu' for other parameters):\n", i + 1, i + 1);
        jabs_message(MSG_INFO, stderr, "  type = %s\n", detector_type_name(det));
        jabs_message(MSG_INFO, stderr, "  theta = %.3lf deg\n", det->theta / C_DEG);
        jabs_message(MSG_INFO, stderr, "  phi = %.3lf deg\n", det->phi / C_DEG);
        if(sim->params->beta_manual) {
            jabs_message(MSG_INFO, stderr, "  beta = %.3lf deg (calculated)\n", sim_exit_angle(sim, det) / C_DEG);
            jabs_message(MSG_INFO, stderr, "  beta = %.3lf deg (manual)\n", det->beta / C_DEG);
        } else {
            jabs_message(MSG_INFO, stderr, "  beta = %.3lf deg\n", sim_exit_angle(sim, det) / C_DEG);
        }
        jabs_message(MSG_INFO, stderr, "  angle from horizontal = %.3lf deg\n", detector_angle(det, 'x')/C_DEG);
        jabs_message(MSG_INFO, stderr, "  angle from vertical = %.3lf deg\n", detector_angle(det, 'y')/C_DEG);
        jabs_message(MSG_INFO, stderr, "  solid angle (given, used) = %.4lf msr\n", i, det->solid/C_MSR);
        if(det->distance > 1.0 * C_MM) {
            jabs_message(MSG_INFO, stderr, "  solid angle (calculated, not used) = %.4lf msr\n", i, detector_solid_angle_calc(det)/C_MSR);
            jabs_message(MSG_INFO, stderr, "  distance = %.3lf mm\n", i, det->distance / C_MM);
            rot_vect v = rot_vect_from_angles(det->theta, det->phi);
            double r = det->distance;
            jabs_message(MSG_INFO, stderr, "  coordinates = (%.3lf, %.3lf, %.3lf) mm\n", v.x * r / C_MM, v.y * r / C_MM, v.z * r / C_MM);
        }
        jabs_message(MSG_INFO, stderr, "  particle solid angle product = %e sr\n", i, sim->fluence * det->solid);
    }
    jabs_message(MSG_INFO, stderr, "n_reactions = %zu\n", sim->n_reactions);
    jabs_message(MSG_INFO, stderr, "fluence = %e (%.5lf p-uC)\n", sim->fluence, sim->fluence*C_E*1.0e6);
    if(sim->channeling_offset != 1.0 || sim->channeling_slope != 0.0) {
        jabs_message(MSG_INFO, stderr, "substrate channeling yield correction offset = %.5lf\n", sim->channeling_offset);
        jabs_message(MSG_INFO, stderr, "substrate channeling yield correction slope = %g / keV (%e)\n", sim->channeling_slope/(1.0/C_KEV), sim->channeling_slope);
    }
}
void sim_workspace_histograms_reset(sim_workspace *ws) {
    for(size_t i = 0; i < ws->n_reactions; i++) {
        sim_reaction *r = ws->reactions[i];
        if(!r)
            continue;
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "Reaction %i:\n", i);
#endif
        gsl_histogram_reset(r->histo);
    }
}


void sim_workspace_histograms_calculate(sim_workspace *ws) {
    for(size_t i = 0; i < ws->n_reactions; i++) {
        sim_reaction *r = ws->reactions[i];
        if(!r)
            continue;
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "Reaction %i:\n", i);
#endif
        bricks_calculate_sigma(ws->det, r->p.isotope, r->bricks, r->last_brick);
        bricks_convolute(r->histo, r->bricks, r->last_brick, ws->fluence * ws->det->solid, ws->params->sigmas_cutoff, ws->params->gaussian_accurate);
    }
}

void sim_workspace_histograms_scale(sim_workspace *ws, double scale) {
    for(size_t i = 0; i < ws->n_reactions; i++) {
        sim_reaction *r = ws->reactions[i];
        if(!r)
            continue;
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "Reaction %i:\n", i);
#endif
        gsl_histogram_scale(r->histo, scale);
    }
}

sim_reaction *sim_reaction_init(sim_workspace *ws, const reaction *r) {
    if(!r) {
        return NULL;
    }
    assert(r->product);
    sim_reaction *sim_r = malloc(sizeof(sim_reaction));
    sim_r->r = r;
    ion *p = &sim_r->p;
    ion_reset(p);
    sim_r->max_depth = 0.0;
    sim_r->i_isotope = ws->sample->n_isotopes; /* Intentionally not valid */

    for(size_t i_isotope = 0; i_isotope < ws->sample->n_isotopes; i_isotope++) {
        if(ws->sample->isotopes[i_isotope] == r->target) {
#ifdef DEBUG
            fprintf(stderr, "Reaction target isotope %s is isotope number %zu in sample.\n", r->target->name, i_isotope);
#endif
            sim_r->i_isotope = i_isotope;
        }
    }
    sim_r->histo = gsl_histogram_alloc(ws->n_channels); /* free'd by sim_workspace_free */
    spectrum_set_calibration(sim_r->histo, ws->det, r->product->Z); /* Setting histogram with Z-specific (or as fallback, default) calibration. */
    gsl_histogram_reset(sim_r->histo);
    sim_r->n_bricks = ws->n_bricks;
    sim_r->bricks = calloc(sim_r->n_bricks, sizeof(brick));
    ion_set_isotope(p, r->product);
    if(p->isotope == ws->ion.isotope) {
        p->nucl_stop = nuclear_stopping_shared_copy(ws->ion.nucl_stop);
    } else {
        p->nucl_stop = nuclear_stopping_new(p->isotope, ws->isotopes);
    }
    sim_reaction_set_cross_section_by_type(sim_r);
    return sim_r;
}

void sim_reaction_free(sim_reaction *sim_r) {
    if(!sim_r) {
        return;
    }
    if(sim_r->histo) {
        gsl_histogram_free(sim_r->histo);
        sim_r->histo = NULL;
    }
    if(sim_r->bricks) {
        free(sim_r->bricks);
        sim_r->bricks = NULL;
    }
    nuclear_stopping_free(sim_r->p.nucl_stop);
    free(sim_r);
}

void sim_reaction_recalculate_internal_variables(sim_reaction *sim_r, double theta, double E_min, double E_max) {
    /* Calculate variables for Rutherford (and Andersen) cross sections. This is done for all reactions, even if they are not RBS or ERD reactions! */
    (void) E_min; /* Energy range could be used to set something (in the future) */
    (void) E_max;
    if(!sim_r || !sim_r->r)
        return;
    const jibal_isotope *incident = sim_r->r->incident;
    const jibal_isotope *target = sim_r->r->target;
    sim_r->E_cm_ratio = target->mass / (incident->mass + target->mass);
    sim_r->mass_ratio = incident->mass / target->mass;
    sim_r->theta = theta;
    if(sim_r->r->Q == 0.0) {
        sim_r->K = reaction_product_energy(sim_r->r, sim_r->theta, 1.0);
    } else {
        sim_r->K = 0.0;
    }
    sim_r->cs_constant = 0.0;
    sim_r->theta_cm = 0.0; /* Will be recalculated, if possible */
    reaction_type type = sim_r->r->type;

    if(!reaction_is_possible(sim_r->r, theta)) {
        sim_r->stop = TRUE;
#ifdef DEBUG
        fprintf(stderr, "Reaction not possible, returning.\n");
#endif
        return;
    }
    if(type == REACTION_RBS) {
        sim_r->theta_cm = sim_r->theta + asin(sim_r->mass_ratio * sin(sim_r->theta));
    }
    if(type == REACTION_RBS_ALT) {
        sim_r->theta_cm = C_PI - (asin(sim_r->mass_ratio * sin(sim_r->theta)) - sim_r->theta);
    }
    if(type == REACTION_RBS || type ==REACTION_RBS_ALT) {
        sim_r->cs_constant = fabs((pow2(sin(sim_r->theta_cm))) / (pow2(sin(sim_r->theta)) * cos(sim_r->theta_cm - sim_r->theta)) *
                                  pow2((incident->Z * C_E * target->Z * C_E) / (4.0 * C_PI * C_EPSILON0)) *
                                  pow4(1.0 / sin(sim_r->theta_cm / 2.0)) * (1.0 / 16.0));
    }
     if(type == REACTION_ERD) { /* ERD */
        sim_r->theta_cm = C_PI - 2.0 * sim_r->theta;
        sim_r->cs_constant = pow2(incident->Z * C_E * target->Z * C_E / (8 * C_PI * C_EPSILON0)) * pow2(1.0 + incident->mass / target->mass) * pow(cos(sim_r->theta), -3.0) * pow2(sim_r->E_cm_ratio);
    }
    if(sim_r->r->cs == JIBAL_CS_ANDERSEN) {
        sim_r->r_VE_factor = 48.73 * C_EV * incident->Z * target->Z * sqrt(pow(incident->Z, 2.0 / 3.0) + pow(target->Z, 2.0 / 3.0)); /* Factors for Andersen correction */
        sim_r->r_VE_factor2 = pow2(0.5 / sin(sim_r->theta_cm / 2.0));
    }
#ifdef DEBUG
    fprintf(stderr, "Reaction recalculated, theta = %g deg, theta_cm = %g deg, K = %g (valid for RBS and ERD). Q = %g MeV.\n", sim_r->theta/C_DEG, sim_r->theta_cm/C_DEG, sim_r->K, sim_r->r->Q / C_MEV);
#endif
}

void sim_reaction_reset_bricks(sim_reaction *sim_r) {
    memset(sim_r->bricks, 0, sizeof(brick) * sim_r->n_bricks);
}

void sim_reaction_set_cross_section_by_type(sim_reaction *sim_r) {
    switch(sim_r->r->type) {
        case REACTION_RBS:
            /* Falls through */
        case REACTION_RBS_ALT:
            /* Falls through */
        case REACTION_ERD:
            sim_r->cross_section = sim_reaction_cross_section_rutherford;
            break;
        case REACTION_FILE:
            sim_r->cross_section = sim_reaction_cross_section_tabulated;
            break;
#ifdef JABS_PLUGINS
        case REACTION_PLUGIN:
            sim_r->cross_section = sim_reaction_cross_section_plugin;
            break;
#endif
        default:
#ifdef DEBUG
            fprintf(stderr, "Unknown reaction type %i, cross_section() function pointer is NULL!\n", sim_r->r->type);
#endif
            sim_r->cross_section = NULL;
    }
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
    if(E > r->E_max || E < r->E_min)
        return 0.0;
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
#ifdef REACTIONS_FALL_BACK
        return sim_reaction_cross_section_rutherford(sim_r, E); /* Fall back quietly to analytical formulae outside tabulated values */
#else
        return 0.0;
#endif
    }
    if(fabs(sim_r->theta - sim_r->r->theta) > 0.01 * C_DEG) {
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "Reaction theta %g deg different from one in the file: %g deg\n", sim_r->theta/C_DEG, sim_r->r->theta/C_DEG);
#endif
#ifdef REACTIONS_FALL_BACK
        return sim_reaction_cross_section_rutherford(sim_r, E);
#else
        return 0.0;
#endif
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

#ifdef JABS_PLUGINS
double sim_reaction_cross_section_plugin(const sim_reaction *sim_r, double E) {
    jabs_plugin_reaction *r = sim_r->r->plugin_r;
    return r->cs(r, sim_r->theta, E);
}
#endif

void sim_sort_reactions(const simulation *sim) {
    qsort(sim->reactions, sim->n_reactions, sizeof(reaction *), &reaction_compare);
}

void sim_reaction_product_energy_and_straggling(sim_reaction *r, const ion *incident) {
    if(r->r->Q == 0.0) {
        r->p.E = incident->E * r->K;
        r->p.S = incident->S * pow2(r->K);
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "Product energy %g keV, eloss straggling %g keV FWHM. Calculated using K = %g\n", r->p.E/C_KEV, C_FWHM * sqrt(r->p.S) / C_KEV, r->K);
#endif
        return;
    }
    r->p.E = reaction_product_energy(r->r, r->theta, incident->E);
    double epsilon = 0.001*C_KEV;
    double deriv = (reaction_product_energy(r->r, r->theta, incident->E+epsilon) - r->p.E)/(epsilon); /* TODO: this derivative could be solved analytically */
    r->p.S = incident->S * pow2(deriv) * incident->E;
#ifdef DEBUG_VERBOSE
    fprintf(stderr, "deriv %g, E_out/E %g, E_out = %g keV, E = %g keV\n", deriv, r->p.E / incident->E, r->p.E/C_KEV, incident->E/C_KEV);
#endif
}

double sim_alpha_angle(const simulation *sim) { /* Returns alpha angle (no sign!) */
    double theta, phi; /* Temporary variables */
    rotate(0.0, 0.0, sim->sample_theta, sim->sample_phi, &theta, &phi); /* Sample in beam system. */
    return theta;
}

double sim_exit_angle(const simulation *sim, const detector *det) { /* Not actually used by simulation, just a convenience function! */
    double theta, phi;
    rotate(sim->sample_theta, sim->sample_phi, det->theta, det->phi, &theta, &phi);
    return C_PI - theta;
}

int sim_do_we_need_erd(const simulation *sim) {
    if(!sim->erd) {
        return FALSE; /* ERD has been (intentionally) disabled */
    }
    if(sim->params->ds) {
        return TRUE; /* In case of DS, ERD becomes possible kind-of possible with any detector geometry */
    }
    int forward_angles = FALSE;
    for(size_t i_det = 0; i_det < sim->n_det; i_det++) {
        const detector *det = sim_det(sim, i_det);
        if(det->theta < C_PI_2) {
            forward_angles = TRUE;
            break;
        }
    }
    return forward_angles; /* If any detector is in forward angle, we might need ERD */
}

void sim_prepare_ion(ion *ion, const simulation *sim, const jibal_isotope *isotopes) {
    ion_reset(ion);
    ion_set_isotope(ion, sim->beam_isotope);
    ion->E = sim->beam_E;
    ion->S = sim->beam_E_broad;
    ion->nucl_stop = nuclear_stopping_new(ion->isotope, isotopes);
}
