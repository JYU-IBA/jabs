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
    sim->params = sim_calc_params_defaults();
    sim->rbs = TRUE;
    sim->erd = TRUE;
    sim->cs_rbs =  jibal->config->cs_rbs;
    sim->cs_erd =  jibal->config->cs_erd;
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
    free(sim);
}

sim_calc_params sim_calc_params_defaults() {
    sim_calc_params p;
    p.stop_step_incident = STOP_STEP_INCIDENT;
    p.stop_step_exiting = STOP_STEP_EXITING;
    p.stop_step_fudge_factor = STOP_STEP_FUDGE_FACTOR;
    p.stop_step_min = STOP_STEP_MIN;
    p.depthsteps_max = 0; /* automatic */
    p.geostragg = FALSE;
    p.beta_manual  = FALSE;
    p.ds = FALSE;
    p.ds_steps_azi = 0;
    p.ds_steps_polar = 0;
    p.rk4 = TRUE;
    p.nucl_stop_accurate = TRUE;
    p.mean_conc_and_energy = FALSE;
    p.cs_stragg_half_n = CS_STRAGG_HALF_N;
    p.cs_n_steps = CS_CONC_STEPS;
    sim_calc_params_update(&p);
    return p;
}

void sim_calc_params_update(sim_calc_params *p) {
    assert(p->cs_n_steps >= 1);
    p->cs_frac = 1.0/(1.0*(p->cs_n_steps+1));
    assert(p.cs_stragg_half_n >= 0);
    p->cs_n_stragg_steps = p->cs_stragg_half_n * 2 + 1;
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
    sim_calc_params_update(p);
}

void sim_calc_params_fast(sim_calc_params *p, int fast) {
    if (fast) {
        p->rk4 = FALSE;
        p->nucl_stop_accurate = FALSE;
        p->mean_conc_and_energy = TRUE;
        p->cs_stragg_half_n = 0; /* Not used if mean_conc_and_energy == TRUE */
        p->cs_n_steps = 0; /* Not used if mean_conc_and_energy == TRUE */
        p->stop_step_fudge_factor *= 1.4;
        p->geostragg = FALSE;
    }
    sim_calc_params_update(p);
}

jibal_cross_section_type sim_cs(const simulation *sim, const reaction_type type) {
    if(type == REACTION_RBS)
        return sim->cs_rbs;
    if(type == REACTION_ERD)
        return sim->cs_erd;
    return JIBAL_CS_NONE;
}

int sim_reactions_add_reaction(simulation *sim, reaction *r) {
    if(!sim || !r)
        return EXIT_FAILURE;
    sim->n_reactions++;
    sim->reactions = realloc(sim->reactions, sim->n_reactions*sizeof(reaction *));
    sim->reactions[sim->n_reactions - 1] = r;
    jabs_message(MSG_INFO, stderr, "Added reaction %zu (%s), %s(%s,%s)%s.\n", sim->n_reactions, reaction_name(r), r->target->name, r->incident->name, r->product->name, r->product_nucleus->name);
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

int sim_reactions_add_auto(simulation *sim, const sample_model *sm, reaction_type type, jibal_cross_section_type cs) { /* Note that sim->ion needs to be set! */
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
    for (size_t i = 0; i < sample->n_isotopes; i++) {
        reaction *r_new = reaction_make(sim->beam_isotope, sample->isotopes[i], type, cs);
        if (!r_new) {
            jabs_message(MSG_ERROR, stderr, "Failed to make an %s reaction with isotope %zu (%s)\n", jibal_cross_section_name(cs), i, sample->isotopes[i]->name);
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
    if(ws->params.stop_step_min == 0.0) { /* Automatic, calculate here. */
        if(det->type == DETECTOR_ENERGY) {
            ws->params.stop_step_min = sqrt(det->calibration->resolution_variance)/2.0;
        } else {
            ws->params.stop_step_min = STOP_STEP_MIN_FALLBACK;
        }
    }
#ifdef DEBUG
    fprintf(stderr, "Minimum stop step %g keV.\n", ws->params.stop_step_min/C_KEV);
#endif
    ws->gsto = jibal->gsto;
    ws->isotopes = jibal->isotopes;
    ws->n_reactions = 0; /* Will be incremented later */

    if(sim->n_reactions == 0) {
        jabs_message(MSG_ERROR, stderr,  "No reactions! Will not initialize workspace if there is nothing to simulate.\n");
        free(ws);
        return NULL;
    }

    if(sim->params.beta_manual && sim->params.ds) {
        jabs_message(MSG_WARNING, stderr,  "Manual exit angle is enabled in addition to dual scattering. This is an unsupported combination. Manual exit angle calculation will be disabled.\n");
        ws->params.beta_manual  = FALSE;
    }

    if(sim->params.beta_manual && sim->params.geostragg) {
        jabs_message(MSG_WARNING, stderr,  "Manual exit angle is enabled in addition to geometric scattering. This is an unsupported combination. Geometric straggling calculation will be disabled.\n");
        ws->params.geostragg = FALSE;
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
    spectrum_set_calibration(ws->histo_sum, ws->det, JIBAL_ANY_Z); /* Calibration (assuming default calibration) can be set now. */
    gsl_histogram_reset(ws->histo_sum); /* This is not necessary, since contents should be set after simulation is over (successfully). */

    size_t n_bricks = 0;

    if(ws->params.depthsteps_max) {
        n_bricks = ws->params.depthsteps_max;
    } else {
        if(det->type == DETECTOR_ENERGY) {
            if(ws->params.stop_step_incident == 0.0) { /* Automatic incident step size */
                n_bricks = (int) ceil(sim->beam_E / (ws->params.stop_step_fudge_factor*sqrt(detector_resolution(ws->det, sim->beam_isotope, sim->beam_E))) + ws->sample->n_ranges);
            } else {
                n_bricks = (int) ceil(sim->beam_E / ws->params.stop_step_incident + ws->sample->n_ranges); /* This is conservative */
            }
        } /* TODO: other detectors? */
        if(n_bricks > 10000) {
            jabs_message(MSG_WARNING, stderr,  "Caution: large number of bricks will be used in the simulation (%zu).\n", n_bricks);
        }
        if(n_bricks > 100000) {
            jabs_message(MSG_WARNING, stderr,  "Number of bricks limited to 100000.\n", n_bricks);
            n_bricks = 100000;
        }
        if(n_bricks < 1000) {
            n_bricks = 1000;
        }
    }
#ifdef DEBUG
    fprintf(stderr, "Number of bricks: %zu\n", n_bricks);
#endif
    ws->reactions = calloc(sim->n_reactions, sizeof (sim_reaction));
    for(size_t i_reaction = 0; i_reaction < sim->n_reactions; i_reaction++) {
        sim_reaction *r = &ws->reactions[ws->n_reactions];
        r->r = sim->reactions[i_reaction];
        if(!r->r) { /* No reaction, this will not be valid */
            continue;
        }
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
        r->histo = gsl_histogram_alloc(ws->n_channels); /* free'd by sim_workspace_free */
        spectrum_set_calibration(r->histo, ws->det, r->r->product->Z); /* Setting histogram with Z-specific (or as fallback, default) calibration. */
        gsl_histogram_reset(r->histo);
        r->n_bricks = n_bricks;
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
        ws->n_reactions++;
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
    double E_max = 0.0;
    for(size_t i_reaction = 0; i_reaction < sim->n_reactions; i_reaction++) {
        double E = reaction_product_energy(sim->reactions[i_reaction], ws->det->theta, sim->beam_E);
        if(E > E_max) {
            E_max = E;
        }
    }
    E_max *= 1.1;
    E_max += detector_resolution(ws->det, sim->beam_isotope, E_max);
#ifdef DEBUG
    fprintf(stderr, "E_max of this simulation is %g keV\n", E_max/C_KEV);
#endif
    size_t i=0;
    while(detector_calibrated(ws->det, JIBAL_ANY_Z, i) < E_max && i <= 1000000) {i++;} /* TODO: this requires changes for ToF spectra. Also TODO: Z-specific calibrations are ignored here */
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
    if(sim->beam_isotope) {
        jabs_message(MSG_INFO, f,  "ion = %s (Z = %i, A = %i, mass %.3lf u)\n", sim->beam_isotope->name, sim->beam_isotope->Z, sim->beam_isotope->A, sim->beam_isotope->mass / C_U);
    } else {
        jabs_message(MSG_INFO, f, "ion = None\n");
    }
    jabs_message(MSG_INFO, f, "E = %.3lf keV\n", sim->beam_E/C_KEV);
    jabs_message(MSG_INFO, f, "E_broad = %.3lf keV FWHM\n", sqrt(sim->beam_E_broad)*C_FWHM/C_KEV);
    jabs_message(MSG_INFO, f, "E_min = %.3lf keV\n", sim->emin/C_KEV);
    jabs_message(MSG_INFO, f, "alpha = %.3lf deg\n", sim_alpha_angle(sim)/C_DEG);
    jabs_message(MSG_INFO, f, "sample tilt (horizontal) = %.3lf deg\n", angle_tilt(sim->sample_theta, sim->sample_phi, 'x')/C_DEG);
    jabs_message(MSG_INFO, f, "sample tilt (vertical) = %.3lf deg\n", angle_tilt(sim->sample_theta, sim->sample_phi, 'y')/C_DEG);
    rot_vect v = rot_vect_from_angles(C_PI - sim->sample_theta, sim->sample_phi); /* By default our sample faces the beam and tilt angles are based on that choice. Pi is there for a reason. */
    jabs_message(MSG_INFO, f, "surf normal unit vector (beam in z direction) = (%.3lf, %.3lf, %.3lf)\n", v.x, v.y, v.z);
    jabs_message(MSG_INFO, f, "aperture = %s ", aperture_name(sim->beam_aperture));
    if(sim->beam_aperture) {
        if(sim->beam_aperture->type == APERTURE_CIRCLE) {
            jabs_message(MSG_INFO, f, "diameter %g mm\n", sim->beam_aperture->diameter/C_MM);
        } else if(sim->beam_aperture->type == APERTURE_RECTANGLE) {
            jabs_message(MSG_INFO, f, "width %g mm height %g mm\n", sim->beam_aperture->width/C_MM, sim->beam_aperture->height/C_MM);
        }
    } else {
        jabs_message(MSG_INFO, f, "\n");
    }
    jabs_message(MSG_INFO, f, "n_detectors = %zu\n", sim->n_det);
    for(size_t i = 0; i < sim->n_det; i++) {
        detector *det = sim->det[i];
        jabs_message(MSG_INFO, f, "DETECTOR %zu (run 'show detector %zu' for other parameters):\n", i + 1, i + 1);
        jabs_message(MSG_INFO, f, "  type = %s\n", detector_type_name(det));
        jabs_message(MSG_INFO, f, "  theta = %.3lf deg\n", det->theta / C_DEG);
        jabs_message(MSG_INFO, f, "  phi = %.3lf deg\n", det->phi / C_DEG);
        if(sim->params.beta_manual) {
            jabs_message(MSG_INFO, f, "  beta = %.3lf deg (calculated)\n", sim_exit_angle(sim, det) / C_DEG);
            jabs_message(MSG_INFO, f, "  beta = %.3lf deg (manual)\n", det->beta / C_DEG);
        } else {
            jabs_message(MSG_INFO, f, "  beta = %.3lf deg\n", sim_exit_angle(sim, det) / C_DEG);
        }
        jabs_message(MSG_INFO, f, "  angle from horizontal = %.3lf deg\n", detector_angle(det, 'x')/C_DEG);
        jabs_message(MSG_INFO, f, "  angle from vertical = %.3lf deg\n", detector_angle(det, 'y')/C_DEG);
        jabs_message(MSG_INFO, f, "  solid angle (given, used) = %.4lf msr\n", i, det->solid/C_MSR);
        if(det->distance > 1.0 * C_MM) {
            jabs_message(MSG_INFO, f, "  solid angle (calculated, not used) = %.4lf msr\n", i, detector_solid_angle_calc(det)/C_MSR);
            jabs_message(MSG_INFO, f, "  distance = %.3lf mm\n", i, det->distance / C_MM);
            rot_vect v = rot_vect_from_angles(det->theta, det->phi);
            double r = det->distance;
            jabs_message(MSG_INFO, f, "  coordinates = (%.3lf, %.3lf, %.3lf) mm\n", v.x * r / C_MM, v.y * r / C_MM, v.z * r / C_MM);
        }
        jabs_message(MSG_INFO, f, "  particle solid angle product = %e sr\n", i, sim->fluence * det->solid);
    }
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
        if(!r)
            continue;
#ifdef DEBUG_VERBOSE
        fprintf(stderr, "Reaction %i:\n", i);
#endif
        brick_int2(r->histo, r->bricks, r->last_brick, ws->det, r->p.isotope, ws->fluence * ws->det->solid);
    }
}

void sim_reaction_recalculate_internal_variables(sim_reaction *sim_r) {
    /* Calculate variables for Rutherford (and Andersen) cross sections. This is done for all reactions, even if they are not RBS or ERD reactions! */
    if(!sim_r || !sim_r->r)
        return;
    const jibal_isotope *incident = sim_r->r->incident;
    const jibal_isotope *target = sim_r->r->target;
    const jibal_isotope *product = sim_r->r->product;
    sim_r->E_cm_ratio = target->mass / (incident->mass + target->mass);
    sim_r->mass_ratio = incident->mass / target->mass;
    sim_r->K = reaction_product_energy(sim_r->r, sim_r->theta, 1.0);
    sim_r->cs_constant = 0.0;
    sim_r->theta_cm = 0.0;
    if(product == incident) { /* RBS */
        if(incident->mass >= target->mass && sim_r->theta > asin(target->mass / incident->mass)) {
            sim_r->stop = TRUE;
            return;
        }
        sim_r->theta_cm = sim_r->theta + asin(sim_r->mass_ratio * sin(sim_r->theta));
        sim_r->cs_constant = (pow2(sin(sim_r->theta_cm))) / (pow2(sin(sim_r->theta)) * cos(sim_r->theta_cm - sim_r->theta)) *
                             pow2((incident->Z * C_E * target->Z * C_E) / (4.0 * C_PI * C_EPSILON0)) *
                             pow4(1.0 / sin(sim_r->theta_cm / 2.0)) * (1.0 / 16.0);
    } else if(product == target) { /* ERD */
        if(sim_r->theta > C_PI/2.0) {
#ifdef DEBUG
            fprintf(stderr, "ERD with %s is not possible (theta %g deg > 90.0 deg)\n", target->name, sim_r->theta/C_DEG);
#endif
            sim_r->stop = TRUE;
            return;
        }
        sim_r->theta_cm = C_PI - 2.0 * sim_r->theta;
        sim_r->cs_constant = pow2(incident->Z * C_E * target->Z * C_E / (8 * C_PI * C_EPSILON0)) * pow2(1.0 + incident->mass / target->mass) * pow(cos(sim_r->theta), -3.0) * pow2(sim_r->E_cm_ratio);
    }
    sim_r->r_VE_factor = 48.73 * C_EV * incident->Z * target->Z * sqrt(pow(incident->Z, 2.0 / 3.0) + pow(target->Z, 2.0 / 3.0)); /* Factors for Andersen correction */
    sim_r->r_VE_factor2 = pow2(0.5 / sin(sim_r->theta_cm / 2.0));
#ifdef DEBUG
    fprintf(stderr, "Reaction recalculated, theta = %g deg, theta_cm = %g deg, K = %g (valid for RBS and ERD). Q = %g MeV.\n", sim_r->theta/C_DEG, sim_r->theta_cm/C_DEG, sim_r->K, sim_r->r->Q / C_MEV);
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
#ifdef REACTIONS_FALL_BACK
    if(E < t[lo].E || E > t[hi].E) {
        return sim_reaction_cross_section_rutherford(sim_r, E); /* Fall back quietly to analytical formulae outside tabulated values */
    }
#endif
    if(fabs(sim_r->theta - sim_r->r->theta) > 0.01 * C_DEG) {
#ifdef DEBUG_VERBOSE
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

void sim_sort_reactions(const simulation *sim) {
    qsort(sim->reactions, sim->n_reactions, sizeof(reaction *), &reaction_compare);
}

void sim_reaction_product_energy_and_straggling(sim_reaction *r, const ion *incident) {
    if(r->r->Q == 0.0) {
        r->p.E = incident->E * r->K;
        r->p.S = incident->S * pow2(r->K);
#ifdef DEBUG
        fprintf(stderr, "Product energy %g keV, eloss straggling %g keV FWHM. Calculated using K = %g\n", r->p.E/C_KEV, C_FWHM * sqrt(r->p.S) / C_KEV, r->K);
#endif
        return;
    }
    r->p.E = reaction_product_energy(r->r, r->theta, incident->E);
    double epsilon = 0.01*C_KEV;
    double deriv = (reaction_product_energy(r->r, r->theta, incident->E+epsilon) - r->p.E)/(epsilon); /* TODO: this derivative could be solved analytically */
    r->p.S = incident->S * pow2(deriv) * incident->E;
#ifdef DEBUG
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
    if(sim->params.ds) {
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
