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
    sim->cs_erd = jabs_reaction_cs_from_jibal_cs(jibal->config->cs_erd);
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

jabs_reaction_cs sim_cs(const simulation *sim, const reaction_type type) {
    if(type == REACTION_RBS || type == REACTION_RBS_ALT)
        return sim->cs_rbs;
    if(type == REACTION_ERD)
        return sim->cs_erd;
    return JABS_CS_NONE;
}

int sim_reactions_add_reaction(simulation *sim, reaction *r) {
    if(!sim || !r)
        return EXIT_FAILURE;
    sim->n_reactions++;
    sim->reactions = realloc(sim->reactions, sim->n_reactions * sizeof(reaction *));
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
                 filename, r->incident->name, r->target->name, r->product->name, r->theta / C_DEG);
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
    if(cs == JABS_CS_NONE) {
        return 0;
    }
    struct sample *sample = sample_from_sample_model(sm);
    if(!sample) {
        return -1;
    }
    for(size_t i = 0; i < sample->n_isotopes; i++) {
        const jibal_isotope *isotope = sample->isotopes[i];
        if(type == REACTION_RBS_ALT) {
            if(sim->beam_isotope->mass <= isotope->mass) { /* alternative solution for RBS exists only if m1 > m2 */
                continue;
            }
        }
        reaction *r_new = reaction_make(sim->beam_isotope, isotope, type, cs);
        if(!r_new) {
            jabs_message(MSG_ERROR, stderr, "Failed to make an %s reaction with %s cross sections for isotope %zu (%s)\n", reaction_type_to_string(type), jabs_reaction_cs_to_string(cs), i,
                         isotope->name);
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

int sim_sanity_check(
        const simulation *sim) { /* This does not guarantee sanity, but it guarantees insanity... Checks should be limited to obvious things. Note that sample and reactions might not be set before this is called. */
    if(!sim) {
        jabs_message(MSG_ERROR, stderr, "No simulation.\n");
        return -1;
    }

    if(!sim->beam_isotope) {
        jabs_message(MSG_ERROR, stderr, "No valid isotope given for the beam.\n");
        return -1;
    }
    if(sim->beam_isotope->Z == 0) {
        jabs_message(MSG_ERROR, stderr, "Neutron beams are not supported.\n");
        return -1;
    }
    if(sim->beam_E > 1000.0 * C_MEV || sim->beam_E < 10 * C_KEV) {
        jabs_message(MSG_ERROR, stderr, "Hmm...? Check your numbers. Your energy is %g J (%g MeV)!\n", sim->beam_E, sim->beam_E / C_MEV);
        return -1;
    }
    if(sim->fluence < 0.0) {
        jabs_message(MSG_ERROR, stderr, "Fluence is negative (%g).\n", sim->fluence);
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

void sim_print(const simulation *sim) {
    if(!sim) {
        return;
    }
    if(sim->beam_isotope) {
        jabs_message(MSG_INFO, stderr, "ion = %s (Z = %i, A = %i, mass %.3lf u)\n", sim->beam_isotope->name, sim->beam_isotope->Z, sim->beam_isotope->A, sim->beam_isotope->mass / C_U);
    } else {
        jabs_message(MSG_INFO, stderr, "ion = None\n");
    }
    jabs_message(MSG_INFO, stderr, "E = %.3lf keV\n", sim->beam_E / C_KEV);
    jabs_message(MSG_INFO, stderr, "E_broad = %.3lf keV FWHM\n", sim->beam_E_broad / C_KEV);
    jabs_message(MSG_INFO, stderr, "E_min = %.3lf keV\n", sim->emin / C_KEV);
    jabs_message(MSG_INFO, stderr, "alpha = %.3lf deg\n", sim_alpha_angle(sim) / C_DEG);
    jabs_message(MSG_INFO, stderr, "sample tilt (horizontal) = %.3lf deg\n", angle_tilt(sim->sample_theta, sim->sample_phi, 'x') / C_DEG);
    jabs_message(MSG_INFO, stderr, "sample tilt (vertical) = %.3lf deg\n", angle_tilt(sim->sample_theta, sim->sample_phi, 'y') / C_DEG);
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
        jabs_message(MSG_INFO, stderr, "  angle from horizontal = %.3lf deg\n", detector_angle(det, 'x') / C_DEG);
        jabs_message(MSG_INFO, stderr, "  angle from vertical = %.3lf deg\n", detector_angle(det, 'y') / C_DEG);
        jabs_message(MSG_INFO, stderr, "  solid angle (given, used) = %.4lf msr\n", i, det->solid / C_MSR);
        if(det->distance > 1.0 * C_MM) {
            jabs_message(MSG_INFO, stderr, "  solid angle (calculated, not used) = %.4lf msr\n", i, detector_solid_angle_calc(det) / C_MSR);
            jabs_message(MSG_INFO, stderr, "  distance = %.3lf mm\n", i, det->distance / C_MM);
            rot_vect v = rot_vect_from_angles(det->theta, det->phi);
            double r = det->distance;
            jabs_message(MSG_INFO, stderr, "  coordinates = (%.3lf, %.3lf, %.3lf) mm\n", v.x * r / C_MM, v.y * r / C_MM, v.z * r / C_MM);
        }
        jabs_message(MSG_INFO, stderr, "  particle solid angle product = %e sr\n", i, sim->fluence * det->solid);
    }
    jabs_message(MSG_INFO, stderr, "n_reactions = %zu\n", sim->n_reactions);
    jabs_message(MSG_INFO, stderr, "fluence = %e (%.5lf p-uC)\n", sim->fluence, sim->fluence * C_E * 1.0e6);
    if(sim->channeling_offset != 1.0 || sim->channeling_slope != 0.0) {
        jabs_message(MSG_INFO, stderr, "substrate channeling yield correction offset = %.5lf\n", sim->channeling_offset);
        jabs_message(MSG_INFO, stderr, "substrate channeling yield correction slope = %g / keV (%e)\n", sim->channeling_slope / (1.0 / C_KEV), sim->channeling_slope);
    }
}

void sim_sort_reactions(const simulation *sim) {
    qsort(sim->reactions, sim->n_reactions, sizeof(reaction *), &reaction_compare);
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
    ion->S = pow2(sim->beam_E_broad / C_FWHM);
    ion->nucl_stop = nuclear_stopping_new(ion->isotope, isotopes);
}
