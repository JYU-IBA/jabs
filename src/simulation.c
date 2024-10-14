/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include "simulation.h"
#include "defaults.h"
#include "rotate.h"
#include "message.h"
#include "geostragg.h"
#include "jabs_debug.h"

simulation *sim_init(jibal *jibal) {
    simulation *sim = calloc(1, sizeof(simulation));
    sim->beam_isotope = jibal_isotope_find(jibal->isotopes, NULL, 2, 4);
    if(!sim->beam_isotope) {
        return NULL;
    }
    sim->beam_aperture = NULL;
    sim->sample_theta = ALPHA_DEFAULT; /* These defaults are for IBM geometry */
    sim->sample_phi = 0.0;
    sim->fluence = FLUENCE_DEFAULT;
    sim->beam_E = ENERGY_DEFAULT;
    sim->emin = E_MIN_DEFAULT;
    sim->beam_E_broad = 0.0;
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
    ion_reset(&sim->ion);
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
    nuclear_stopping_free(sim->ion.nucl_stop);
    ion_gsto_free(sim->ion.ion_gsto);
    free(sim);
}

jabs_reaction_cs sim_cs(const simulation *sim, const reaction_type type) {
    if(type == REACTION_RBS || type == REACTION_RBS_ALT) {
        return sim->cs_rbs;
    }
    if(type == REACTION_ERD) {
        return sim->cs_erd;
    }
    return JABS_CS_NONE;
}

reaction *sim_reaction_make_from_argv(const jibal *jibal, const simulation *sim, int *argc, char * const **argv) {
    reaction *r = reaction_make_from_argv(jibal, sim->beam_isotope, argc, argv);
    if(!r) {
        return NULL;
    }
    if(r->cs == JABS_CS_NONE) {
        jabs_reaction_cs cs = sim_cs(sim, r->type);
        r->cs = cs;
        jabs_message(MSG_DEBUG, "Reaction cross section not given, assuming default for %s: %s.\n", reaction_type_to_string(r->type), jabs_reaction_cs_to_string(cs));
    }
    return r;
}

int sim_reactions_add_reaction(simulation *sim, reaction *r, int silent) {
    if(!sim || !r)
        return EXIT_FAILURE;
    if(r->product == NULL) {
        jabs_message(MSG_ERROR, "Reaction product is not defined, can not add reaction.\n");
        return EXIT_FAILURE;
    }
    sim->n_reactions++;
    sim->reactions = realloc(sim->reactions, sim->n_reactions * sizeof(reaction *));
    sim->reactions[sim->n_reactions - 1] = r;
    if(!silent) {
        jabs_message(MSG_INFO, "Added reaction %zu (%s), E = [%g keV, %g keV], Q = %g keV\n",
                     sim->n_reactions, reaction_name(r), r->E_min / C_KEV, r->E_max / C_KEV, r->Q / C_KEV);
    }
    return EXIT_SUCCESS;
}

int sim_reactions_remove_reaction(simulation *sim, size_t i) {
    if(i >= sim->n_reactions) {
        jabs_message(MSG_ERROR, "There are %zu reactions, can't remove reaction %zu\n", sim->n_reactions, i + 1);
        return EXIT_FAILURE;
    }
    reaction_free(sim->reactions[i]);
    sim->reactions[i] = NULL;
    sim_sort_reactions(sim);
    sim->n_reactions--;
    return EXIT_SUCCESS;
}

int sim_reactions_add_auto(simulation *sim, const sample_model *sm, reaction_type type, jabs_reaction_cs cs, int silent) { /* Note that sim->ion needs to be set! */
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
            jabs_message(MSG_ERROR, "Failed to make an %s reaction with %s cross sections for isotope %zu (%s)\n",
                         reaction_type_to_string(type), jabs_reaction_cs_to_string(cs), i, isotope->name);
            continue;
        }
        sim_reactions_add_reaction(sim, r_new, silent);
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
        jabs_message(MSG_ERROR, "No simulation.\n");
        return -1;
    }

    if(!sim->beam_isotope) {
        jabs_message(MSG_ERROR, "No valid isotope given for the beam.\n");
        return -1;
    }
    if(sim->beam_isotope->Z == 0) {
        jabs_message(MSG_ERROR, "Neutron beams are not supported.\n");
        return -1;
    }
    if(sim->beam_E >= E_MAX || sim->beam_E <= E_MIN || sim->beam_E <= sim->emin) {
        jabs_message(MSG_ERROR, "Beam energy is %g MeV, it should be < %g MeV and > %g keV and > simulation minimum energy, which is currently set to %g keV.\n", sim->beam_E / C_MEV, E_MAX / C_MEV, E_MIN / C_KEV, sim->emin / C_KEV);
        return -1;
    }
    if(sim->fluence < 0.0) {
        jabs_message(MSG_ERROR, "Fluence is negative (%g).\n", sim->fluence);
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
    if(det->name == NULL) {
        if(asprintf(&(det->name), "Detector %zu", sim->n_det) < 0) {
            det->name = NULL;
        }
    }
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

void sim_print(const simulation *sim, jabs_msg_level msg_level) {
    if(!sim) {
        return;
    }
    if(sim->beam_isotope) {
        jabs_message(msg_level, "ion = %s (Z = %i, A = %i, mass %.3lf u)\n", sim->beam_isotope->name, sim->beam_isotope->Z, sim->beam_isotope->A, sim->beam_isotope->mass / C_U);
    } else {
        jabs_message(msg_level, "ion = None\n");
    }
    jabs_message(msg_level, "E = %.3lf keV\n", sim->beam_E / C_KEV);
    jabs_message(msg_level, "E_broad = %.3lf keV FWHM\n", sim->beam_E_broad / C_KEV);
    jabs_message(msg_level, "E_min = %.3lf keV\n", sim->emin / C_KEV);
    jabs_message(msg_level, "alpha = %.3lf deg\n", sim_alpha_angle(sim) / C_DEG);
    jabs_message(msg_level, "sample tilt (horizontal) = %.3lf deg\n", angle_tilt(sim->sample_theta, sim->sample_phi, 'x') / C_DEG);
    jabs_message(msg_level, "sample tilt (vertical) = %.3lf deg\n", angle_tilt(sim->sample_theta, sim->sample_phi, 'y') / C_DEG);
    rot_vect v = rot_vect_from_angles(C_PI - sim->sample_theta, sim->sample_phi); /* By default our sample faces the beam and tilt angles are based on that choice. Pi is there for a reason. */
    jabs_message(msg_level, "surf normal unit vector (beam in z direction) = (%.3lf, %.3lf, %.3lf)\n", v.x, v.y, v.z);
    char *aperture_str = aperture_to_string(sim->beam_aperture);
    jabs_message(msg_level, "aperture = %s\n", aperture_str);
    free(aperture_str);
    jabs_message(msg_level, "n_detectors = %zu\n", sim->n_det);
    for(size_t i = 0; i < sim->n_det; i++) {
        detector *det = sim->det[i];
        jabs_message(msg_level, "DETECTOR %zu (run 'show detector %zu' for other parameters):\n", i + 1, i + 1);
        jabs_message(msg_level, "  type = %s\n", detector_type_name(det));
        jabs_message(msg_level, "  theta = %.3lf deg\n", det->theta / C_DEG);
        jabs_message(msg_level, "  phi = %.3lf deg\n", det->phi / C_DEG);
        if(sim->params->beta_manual) {
            jabs_message(msg_level, "  beta = %.3lf deg (calculated)\n", sim_exit_angle(sim, det) / C_DEG);
            jabs_message(msg_level, "  beta = %.3lf deg (manual)\n", det->beta / C_DEG);
        } else {
            jabs_message(msg_level, "  beta = %.3lf deg\n", sim_exit_angle(sim, det) / C_DEG);
        }
        jabs_message(msg_level, "  angle from horizontal = %.3lf deg\n", detector_angle(det, 'x') / C_DEG);
        jabs_message(msg_level, "  angle from vertical = %.3lf deg\n", detector_angle(det, 'y') / C_DEG);
        jabs_message(msg_level, "  solid angle (given, used) = %.4lf msr\n", det->solid / C_MSR);
        if(det->distance > 1.0 * C_MM) {
            jabs_message(msg_level, "  solid angle (calculated, not used) = %.4lf msr\n", detector_solid_angle_calc(det) / C_MSR);
            jabs_message(msg_level, "  distance = %.3lf mm\n", det->distance / C_MM);
            rot_vect v = rot_vect_from_angles(det->theta, det->phi);
            double r = det->distance;
            jabs_message(msg_level, "  coordinates = (%.3lf, %.3lf, %.3lf) mm\n", v.x * r / C_MM, v.y * r / C_MM, v.z * r / C_MM);
        }
        jabs_message(msg_level, "  particle solid angle product = %e sr\n", sim->fluence * det->solid);
    }
    jabs_message(msg_level, "n_reactions = %zu\n", sim->n_reactions);
    jabs_message(msg_level, "fluence = %e (%.5lf p-uC)\n", sim->fluence, sim->fluence * C_E * 1.0e6);
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
    return exit_angle(sim->sample_theta, sim->sample_phi, det->theta, det->phi);
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

int sim_prepare_ion(ion *ion, const simulation *sim, const jibal_isotope *isotopes, const jibal_gsto *gsto) {
    nuclear_stopping_free(ion->nucl_stop);
    ion_gsto_free(ion->ion_gsto);
    ion_reset(ion);
    ion_set_isotope(ion, sim->beam_isotope);
    ion->nucl_stop = nuclear_stopping_new(sim->beam_isotope, isotopes);
    ion->ion_gsto = ion_gsto_new(sim->beam_isotope, gsto);
    if(!ion->nucl_stop || !ion->ion_gsto) {
        DEBUGMSG("Preparing ion failed, nucl_stop = %p, ion_gsto = %p.", (void *)ion->nucl_stop, (void *)ion->ion_gsto);
        nuclear_stopping_free(ion->nucl_stop);
        ion_gsto_free(ion->ion_gsto);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int sim_prepare_reactions(const simulation *sim, const jibal_isotope *isotopes, const jibal_gsto *gsto) {
    for(size_t i = 0; i < sim->n_reactions; i++) {
        reaction *r = sim->reactions[i];
        nuclear_stopping_free(r->nucl_stop);
        ion_gsto_free(r->ion_gsto);
        if(r->product == sim->beam_isotope) {
            r->nucl_stop = nuclear_stopping_shared_copy(sim->ion.nucl_stop);
            r->ion_gsto = ion_gsto_shared(sim->ion.ion_gsto);
        } else {
            r->nucl_stop = nuclear_stopping_new(r->product, isotopes);
            r->ion_gsto = ion_gsto_new(r->product, gsto);
        }
        if(!r->nucl_stop || !r->ion_gsto) {
            jabs_message(MSG_ERROR, "Error in preparing reaction %s: problems assigning stopping data.", r->name);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
