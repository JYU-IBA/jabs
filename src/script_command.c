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
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
#include <direct.h> // getcwd
#else
#include <unistd.h>
#endif
#include <jibal_generic.h>
#include <jibal_defaults.h>
#include <gsl/gsl_version.h>
#include "jabs_debug.h"
#include "message.h"
#include "sample.h"
#include "spectrum.h"
#include "simulation.h"
#include "fit.h"
#include "options.h"
#include "git.h"
#include "jabs.h"
#include "generic.h"
#include "script_session.h"
#include "script_command.h"
#include "plugin.h"
#include "idf2jbs.h"
#include "simulation2idf.h"


int script_prepare_sim_or_fit(script_session *s) {
    fit_data *fit = s->fit;
    if(!fit->sm) {
        jabs_message(MSG_ERROR, "No sample has been defined!\n");
        return -1;
    }
    if(!fit->sim->beam_isotope) {
        jabs_message(MSG_ERROR, "No ion has been defined!\n");
        return -1;
    }
    if(!fit->sim->det || fit->sim->n_det == 0) {
        jabs_message(MSG_ERROR, "No detector has been defined!\n");
        return -1;
    }
    if(sim_sanity_check(fit->sim)) {
        jabs_message(MSG_ERROR, "Simulation failed sanity check.\n");
        return -1;
    }
    sample_free(fit->sim->sample);
    fit->sim->sample = sample_from_sample_model(fit->sm);
    if(!fit->sim->sample) {
        jabs_message(MSG_ERROR, "Could not make a sample based on model description. This should never happen.\n");
        return -1;
    }
    for(size_t i = 0; i < fit->sim->n_reactions; i++) {
        const reaction *r = fit->sim->reactions[i];
        if(r->incident != fit->sim->beam_isotope) {
            jabs_message(MSG_INFO, "Reaction %i is for an incident %s, but this simulation should be for %s. Try \"reset reactions\".\n", i + 1, r->incident->name, fit->sim->beam_isotope);
            return EXIT_FAILURE;
        }
    }

    if(fit->sim->n_reactions == 0) {
        jabs_message(MSG_WARNING, "No reactions defined. Please be aware there are commands called \"reset reactions\" and \"add reactions\".\n");
        if(fit->sim->rbs) {
            jabs_message(MSG_INFO, "Adding RBS reactions.\n");
            sim_reactions_add_auto(fit->sim, fit->sm, REACTION_RBS, sim_cs(fit->sim, REACTION_RBS), TRUE);
            sim_reactions_add_auto(fit->sim, fit->sm, REACTION_RBS_ALT, sim_cs(fit->sim, REACTION_RBS_ALT), TRUE);
            /* TODO: loop over all detectors and add reactions that are possible (one reaction for all detectors) */
        }
        if(sim_do_we_need_erd(fit->sim)) {
            jabs_message(MSG_INFO, "Adding ERDA reactions.\n");
            sim_reactions_add_auto(fit->sim, fit->sm, REACTION_ERD, sim_cs(fit->sim, REACTION_ERD), TRUE);
        }
    }
    if(fit->sim->n_reactions == 0) {
        jabs_message(MSG_ERROR, "No reactions. Nothing to do.\n");
        return EXIT_FAILURE;
    }
    sim_sort_reactions(fit->sim);
    jabs_message(MSG_INFO, "Sample model:\n");
    sample_model_print(NULL, fit->sm, MSG_INFO);
    jabs_message(MSG_VERBOSE, "Simplified sample model for simulation:\n");
    sample_print(fit->sim->sample, TRUE, MSG_VERBOSE);

    if(assign_stopping(fit->jibal->gsto, fit->sim)) {
        jabs_message(MSG_ERROR, "Could not assign stopping or straggling. Failure. Provide more data, check that JIBAL Z2_max is sufficiently large (currently %i) or disable unwanted reactions (e.g. ERD).\n",
                     s->jibal->config->Z_max);
        return -1;
    }
    script_show_stopping(s, 0, NULL);
#ifdef DEBUG
    jibal_gsto_print_files(fit->jibal->gsto, TRUE);
#endif
    jabs_message(MSG_VERBOSE, "Loading stopping data.\n");
    jibal_gsto_load_all(fit->jibal->gsto);
    DEBUGSTR("Updating calculation params before sim/fit");
    sim_calc_params_update(fit->sim->params);
    jabs_message(MSG_VERBOSE, "Simulation parameters:\n");
    sim_print(fit->sim, MSG_VERBOSE);

    if(sim_prepare_ion(&fit->sim->ion, fit->sim, fit->jibal->isotopes, fit->jibal->gsto)) {
        jabs_message(MSG_ERROR, "Preparing incident ion failed. Possibly issue with either nuclear stopping or GSTO stopping.\n");
        return EXIT_FAILURE;
    }
    if(fit->sim->ion.ion_gsto->emax < fit->sim->beam_E) {
        jabs_message(MSG_ERROR, "Maximum energy for incident %s is %g keV for some target element(s) based on stopping data (GSTO), which is more than beam energy of %g keV.\n",
                     fit->sim->ion.isotope->name, fit->sim->ion.ion_gsto->emax / C_KEV, fit->sim->beam_E / C_KEV);
        return EXIT_FAILURE;
    }
    if(fit->sim->ion.ion_gsto->emin > fit->sim->emin) {
        jabs_message(MSG_WARNING, "Minimum energy for incident %s is %g keV for some target element(s) based on stopping data (GSTO), which is more than current workspace minimum energy of %g keV.\n",
                     fit->sim->ion.isotope->name, fit->sim->ion.ion_gsto->emin / C_KEV, fit->sim->emin / C_KEV);
    }
    if(sim_prepare_reactions(fit->sim, fit->jibal->isotopes, fit->jibal->gsto)) {
        return EXIT_FAILURE;
    }
    reactions_print(fit->sim->reactions, fit->sim->n_reactions);
    fit_data_spectra_alloc(fit);
    s->start = jabs_clock();
    return 0;
}

int script_finish_sim_or_fit(script_session *s) {
    s->end = jabs_clock();
    double time = s->end - s->start;
    if(time > 1.0) {
        jabs_message(MSG_IMPORTANT, "\n...finished! Total time: %.3lf s.\n", time);
    } else {
        jabs_message(MSG_IMPORTANT, "\n...finished! Total time: %.3lf ms.\n", time * 1000.0);
    }
    fit_data_fdd_free(s->fit);
#ifdef CLEAR_GSTO_ASSIGNMENTS_WHEN_FINISHED
    jibal_gsto_assign_clear_all(s->fit->jibal->gsto); /* Is it necessary? No. Here? No. Does it clear old stuff? Yes. */
#endif
    return 0;
}


void script_command_not_found(const char *cmd, const script_command *c_parent) {
    if(c_parent) {
        if(cmd) {
            jabs_message(MSG_ERROR, "Sub-command \"%s\" is invalid!\n\n", cmd);
        } else {
            jabs_message(MSG_ERROR, "Not enough arguments!\n");
        }
        if(c_parent->subcommands) {
            size_t matches = script_command_print_possible_matches_if_ambiguous(c_parent->subcommands, cmd);
            if(matches == 0) {
                jabs_message(MSG_INFO, "\nFollowing subcommands of \"%s\" are recognized:\n", c_parent->name);
                script_commands_print(c_parent->subcommands);
            }
        }
    } else {
        if(cmd) {
            jabs_message(MSG_ERROR, "Invalid command or argument: \"%s\".\n", cmd);
        } else {
            jabs_message(MSG_ERROR, "What?\n", cmd);
        }
    }
}

script_command_status script_simulate(script_session *s, int argc, char *const *argv) {
    const int argc_orig = argc;
    (void) argv;
    struct fit_data *fit = s->fit;
    if(argc > 1) {
        /* TODO? */
    }
    if(script_prepare_sim_or_fit(s)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    sim_calc_params_print(fit->sim->params, MSG_VERBOSE);
    jabs_message(MSG_IMPORTANT, "Simulation begins...\n");
    for(size_t i_det = 0; i_det < fit->sim->n_det; i_det++) {
        detector *det = fit->sim->det[i_det];
        detector_update(det);
        sim_workspace *ws = sim_workspace_init(s->jibal, fit->sim, det);
        if(simulate_with_ds(ws)) {
            jabs_message(MSG_ERROR, "Simulation failed.\n");
            sim_workspace_free(ws);
            return SCRIPT_COMMAND_FAILURE;
        }
#ifdef DEBUG
        char *bricks_filename;
        asprintf(&bricks_filename, "bricks_%zu.dat", i_det + 1);
        if(bricks_filename) {
            sim_workspace_print_bricks(ws, bricks_filename);
            free(bricks_filename);
        }
#endif
        fit_data_spectra_copy_to_spectra_from_ws(&fit->spectra[i_det], det, s->fit->exp[i_det], ws);
        sim_workspace_free(ws);
    }
    script_finish_sim_or_fit(s);
    return argc_orig - argc;
}

script_command_status script_fit(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit_data = s->fit;
    const char *fit_usage = "Usage: fit [fitvar1,fitvar2,...]\nSee 'show fit variables' for a list of possible fit variables.\n";
    if(argc != 1) {
        jabs_message(MSG_ERROR, fit_usage);
        return SCRIPT_COMMAND_FAILURE;
    }

    fit_params *p_all = fit_params_all(fit_data);
    if(fit_params_enable_using_string(p_all, argv[0])) {
        jabs_message(MSG_ERROR, "Error in adding fit parameters.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(fit_params_update(p_all)) {
        jabs_message(MSG_ERROR, "Error in checking fit parameters.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    fit_params_print(p_all, TRUE, NULL, MSG_VERBOSE);
    fit_params_free(fit_data->fit_params);
    fit_data->fit_params = NULL;
    fit_data->fit_params = p_all;

    if(fit_data->fit_params->n_active == 0) {
        jabs_message(MSG_ERROR, fit_usage);
        return SCRIPT_COMMAND_FAILURE;
    }

    jabs_message(MSG_INFO, "%zu fit parameters possible, %zu active.\n", fit_data->fit_params->n, fit_data->fit_params->n_active);
    if(!fit_data->exp) { /* TODO: not enough to check this */
        jabs_message(MSG_ERROR, "No experimental spectrum set.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(fit_data->n_fit_ranges == 0) {
        jabs_message(MSG_ERROR, "No fit range(s) given. Use 'add fit range'.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(script_prepare_sim_or_fit(s)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    fit_data->fit_iter_callback = s->fit_iter_callback;
    if(fit(fit_data) < 0) {
        return SCRIPT_COMMAND_FAILURE;
    }
    jabs_message(MSG_VERBOSE, "\nFinal profile:\n");
    sample_print(fit_data->sim->sample, FALSE, MSG_VERBOSE);
    jabs_message(MSG_VERBOSE, "\nFinal layer thicknesses:\n");
    sample_print_thicknesses(NULL, fit_data->sim->sample, MSG_VERBOSE);
    jabs_message(MSG_VERBOSE, "\nFinal sample model:\n");
    sample_model_print(NULL, fit_data->sm, MSG_VERBOSE);
    jabs_message(MSG_INFO, "\n");
    fit_stats_print(&fit_data->stats, MSG_INFO);
    script_finish_sim_or_fit(s);
    return 1;
}

script_command_status script_save_simulation(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    const int argc_orig = argc;
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Usage: save simulation <file>\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(simulation2idf(fit, argv[0])) {
        jabs_message(MSG_ERROR, "Could not save simulation to file %s.", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    argc -= 1;
    return argc_orig - argc;
}

script_command_status script_save_spectra(script_session *s, int argc, char *const *argv) {
    size_t i_det = 0;
    struct fit_data *fit = s->fit;
    const int argc_orig = argc;
    if(script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det) || argc < 1) {
        jabs_message(MSG_ERROR, "Usage: save spectra {<detector>} file\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Not enough arguments for save spectra.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(sim_workspace_print_spectra(fit->spectra, argv[0])) {
        jabs_message(MSG_ERROR,
                     "Could not save spectra of detector %zu to file \"%s\"! There should be %zu detector(s).\n",
                     i_det + 1, argv[0], fit->sim->n_det);
        return SCRIPT_COMMAND_FAILURE;
    }
    argc -= 1;
    return argc_orig - argc;
}


script_command_status script_save_sample(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit_data = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Usage: save sample <file>\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(!fit_data->sm) {
        jabs_message(MSG_ERROR, "No sample set.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(sample_model_print(argv[0], fit_data->sm, MSG_INFO)) {
        jabs_message(MSG_ERROR, "Could not write sample to file \"%s\".\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    return 1;
}

script_command_status script_save_calibrations(script_session *s, int argc, char *const *argv) {
    const int argc_orig = argc;
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Usage: save calibrations <file>\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    FILE *f = fopen_file_or_stream(argv[0], "w");
    if(!f) {
        jabs_message(MSG_ERROR, "Can not open file \"%s\" for writing.\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    argc--;
    for(size_t i = 0; i < s->fit->sim->n_det; i++) {
        const detector *det = sim_det(s->fit->sim, i);
        char *calib_str = calibration_to_string(det->calibration);
        char *reso_str = detector_resolution_to_string(det, JIBAL_ANY_Z);
        jabs_message_printf(MSG_INFO, f, "set detector %zu calibration %s resolution %s\n", i + 1, calib_str, reso_str);
        free(calib_str);
        free(reso_str);
    }
    fclose_file_or_stream(f);
    return argc_orig - argc;
}

script_command_status script_remove_reaction(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    static const char *remove_reaction_usage = "Usage: remove reaction {<number | <type> <target isotope>}\n";
    if(argc < 1) {
        jabs_message(MSG_ERROR, remove_reaction_usage);
        return SCRIPT_COMMAND_FAILURE;
    }
    char *end;
    size_t i_reaction = strtoull(argv[0], &end, 10);
    if(*end == '\0') {
        if(sim_reactions_remove_reaction(fit->sim, i_reaction - 1)) {
            return SCRIPT_COMMAND_FAILURE;
        } else {
            return 1; /* Number of consumed arguments */
        }
    }
    if(argc < 2) {
        jabs_message(MSG_ERROR, remove_reaction_usage);
        return SCRIPT_COMMAND_FAILURE;
    }
    reaction_type type = reaction_type_from_string(argv[0]);
    const jibal_isotope *target = jibal_isotope_find(fit->jibal->isotopes, argv[1], 0, 0);
    if(type == REACTION_NONE) {
        jabs_message(MSG_ERROR, "This is not a valid reaction type: \"%s\".\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    if(!target) {
        jabs_message(MSG_ERROR, "This is not a valid isotope: \"%s\".\n", argv[1]);
        return SCRIPT_COMMAND_FAILURE;
    }
    for(size_t i = 0; i < fit->sim->n_reactions; i++) {
        if(fit->sim->reactions[i]->type == type && fit->sim->reactions[i]->target == target) {
            if(sim_reactions_remove_reaction(fit->sim, i)) {
                return SCRIPT_COMMAND_FAILURE;
            } else {
                return 2; /* Number of consumed arguments. */
            }
        }
    }
    jabs_message(MSG_ERROR, "No matching reaction found!\n");
    return SCRIPT_COMMAND_FAILURE;
}

script_command_status script_roi(script_session *s, int argc, char *const *argv) {
    const struct fit_data *fit = s->fit;
    const int argc_orig = argc;
    size_t i_det = 0;

    if(script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det)) {
        return SCRIPT_COMMAND_FAILURE;
    }

    if(argc < 1) {
        jabs_message(MSG_ERROR, "Usage: roi <range> {<range> <range> ...}\nExample: roi [400:900] [980:1200]\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    while(argc > 0) {
        roi r = {.i_det = i_det};
        if(fit_set_roi_from_string(&r, argv[0])) {
            return SCRIPT_COMMAND_FAILURE;
        }
        fit_data_roi_print(s->fit, &r);
        argc--;
        argv++;
    }
    return argc_orig - argc;
}

script_command_status script_exit(script_session *s, int argc, char *const *argv) {
    (void) s;
    (void) argc;
    (void) argv;
    return SCRIPT_COMMAND_EXIT;
}

const script_command *script_command_find(const script_command *commands, const char *cmd_string) {
    if(!cmd_string) {
        return NULL;
    }
    if(*cmd_string == '\0' || *cmd_string == '#') {
        return NULL;
    }
    int found = 0;
    const script_command *c_found = NULL;
    for(const script_command *c = commands; c; c = c->next) {
        if(strncmp(c->name, cmd_string, strlen(cmd_string)) == 0) {
            found++;
            c_found = c;
            DEBUGMSG("Candidate for \"%s\": \"%s\".", cmd_string, c->name);
            if(strlen(cmd_string) == strlen(c->name)) { /* Exact match, can not be ambiguous */
                found = 1;
                break;
            }
        }
    }
    if(found == 1) {
        return c_found;
    }
    return NULL;
}

size_t script_command_print_possible_matches_if_ambiguous(const script_command *commands, const char *cmd_string) {
    size_t found = 0;
    if(!cmd_string || !commands)
        return 0;
    for(const script_command *c = commands; c; c = c->next) {
        if(strncmp(c->name, cmd_string, strlen(cmd_string)) == 0) {
            found++;
            if(strlen(cmd_string) == strlen(c->name)) { /* Exact match, can not be ambiguous */
                found = 1;
                break;
            }
        }
    }
    if(found > 1) {
        jabs_message(MSG_ERROR, "\"%s\" is ambiguous (%i matches):", cmd_string, found);
        for(const script_command *c = commands; c; c = c->next) {
            if(strncmp(c->name, cmd_string, strlen(cmd_string)) == 0) {
                jabs_message(MSG_ERROR, " %s", c->name);
            }
        }
        jabs_message(MSG_ERROR, "\n");
    }
    return found;
}

script_command_status script_set_var(struct script_session *s, jibal_config_var *var, int argc, char *const *argv) {
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Not enough arguments to set variable.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(var->type == JIBAL_CONFIG_VAR_UNIT) { /* Provides error checking, which JIBAL doesn't */
        double out;
        if(jabs_unit_convert(s->jibal->units, var->unit_type, argv[0], &out) < 0) {
            return EXIT_FAILURE;
        }
        *((double *)var->variable) = out;
        return 1;
    }
    if(var->type == JIBAL_CONFIG_VAR_OPTION) { /* Provides slightly different functionality than JIBAL (partial matches are ok) */
        int value;
        size_t l = strlen(argv[0]);
        int found = 0;
        for(const jibal_option *o = var->option_list; o->s; o++) {
            if(strncmp(o->s, argv[0], l) == 0) { /* Partial matches are suitable candidates */
                value = o->val;
                found++;
            }
        }
        if(found == 1) {
            *((int *)var->variable) = value; /* Exactly one hit, assign. */
            return 1;
        }
        if(found > 1) {
            jabs_message(MSG_ERROR, "Value \"%s\" is ambiguous (%i matches).\n", argv[0], found);
        } else if(found == 0) {
            jabs_message(MSG_ERROR, "Value \"%s\" doesn't match with any known option.\n", argv[0]);
        }
        jabs_message(MSG_ERROR, "Valid options are:\n");
        for(const jibal_option *o = var->option_list; o->s; o++) {
            jabs_message(MSG_ERROR, " %s", o->s);
        }
        jabs_message(MSG_ERROR, "\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    jibal_config_var_set(s->jibal->units, var, argv[0], NULL);
    return 1; /* Number of arguments */
}

script_command_status script_set_detector_val(struct script_session *s, int val, int argc, char *const *argv) {
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Not enough arguments to set detector variable.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    DEBUGMSG("Active detector is %zu.", s->i_det_active);
    detector *det = sim_det(s->fit->sim, s->i_det_active);
    if(!det) {
        jabs_message(MSG_ERROR, "Detector %zu does not exist.\n", s->i_det_active);
        return SCRIPT_COMMAND_FAILURE;
    }
    DEBUGMSG("Some kind of value to be converted: \"%s\"", argv[0]);
    double *value_double = NULL;
    size_t *value_size = NULL;
    char unit_type = JIBAL_UNIT_TYPE_ANY;
    switch(val) {
        case 'b': /* beta */
            value_double = &(det->beta);
            unit_type = JIBAL_UNIT_TYPE_ANGLE;
            break;
        case 'c': /* column */
            value_size = &(det->column);
            break;
        case 'h': /* channels */
            value_size = &(det->channels);
            break;
        case 'C': /* compress */
            value_size = &(det->compress);
            break;
        case 'd': /* distance */
            value_double = &(det->distance);
            unit_type = JIBAL_UNIT_TYPE_DISTANCE;
            break;
        case 't': /* type */
            det->type = jibal_option_get_value(detector_option, argv[0]);
            if(det->type == 0) {
                jabs_message(MSG_ERROR, "Detector type \"%s\" is none or unknown.\n", argv[0]);
                return SCRIPT_COMMAND_FAILURE;
            }
            break;
        case 'S': /* slope, this is for backwards compatibility (and ease of use with linear calibration) */
            value_double = calibration_get_param_ref(det->calibration, CALIBRATION_PARAM_SLOPE);
            unit_type = JIBAL_UNIT_TYPE_ANY; /* Could be time or energy */
            break;
        case 'O': /* offset, this is for backwards compatibility */
            value_double = calibration_get_param_ref(det->calibration, CALIBRATION_PARAM_OFFSET);
            unit_type = JIBAL_UNIT_TYPE_ANY;
            break;
        case 'r': /* resolution */
            value_double = calibration_get_param_ref(det->calibration, CALIBRATION_PARAM_RESOLUTION);
            unit_type = JIBAL_UNIT_TYPE_ANY;
            break;
        case 's': /* solid */
            value_double = &(det->solid);
            unit_type = JIBAL_UNIT_TYPE_SOLID_ANGLE;
            break;
        case 'T': /* theta */
            value_double = &(det->theta);
            unit_type = JIBAL_UNIT_TYPE_ANGLE;
            break;
        case 'l': /* length */
            value_double = &(det->length);
            unit_type = JIBAL_UNIT_TYPE_DISTANCE;
            break;
        case 'p': /* phi */
            value_double = &(det->phi);
            unit_type = JIBAL_UNIT_TYPE_ANGLE;
            break;
        default:
            jabs_message(MSG_ERROR, "Unhandled value %i in script_set_detector_val. Report to developer.\n", val);
            return SCRIPT_COMMAND_FAILURE;
            break;
    }
    if(value_double) {
        if(jabs_unit_convert(s->jibal->units, unit_type, argv[0], value_double) < 0) {
            return SCRIPT_COMMAND_FAILURE;
        }
    } else if(value_size) {
        char *end;
        size_t val_converted_ull = strtoull(argv[0], &end, 10);
        if(*end == '\0') {
            *value_size = val_converted_ull;
        } else {
            jabs_message(MSG_ERROR, "Conversion of \"%s\" to unsigned integer failed.\n", argv[0]);
        }
    }
    return 1; /* Number of arguments */
}

script_command_status script_set_detector_calibration_val(struct script_session *s, int val, int argc, char *const *argv) {
    (void) argv;
    DEBUGMSG("Active detector is %zu.", s->i_det_active);
    detector *det = sim_det(s->fit->sim, s->i_det_active);
    if(!det) {
        jabs_message(MSG_ERROR, "Detector %zu does not exist.\n", s->i_det_active);
        return SCRIPT_COMMAND_FAILURE;
    }
    calibration *c;
    switch(val) { /* Handle cases where we don't expect (consume) arguments */
        case 'L': /* linear */
            c = calibration_init_linear();
            calibration_copy_params(c, detector_get_calibration(det, s->Z_active)); /* Copy parameters from old calibration, as much as possible */
            if(detector_set_calibration_Z(s->jibal->config, det, c, s->Z_active)) {
                jabs_message(MSG_ERROR, "Could not set linear calibration (element = %s).\n",
                             jibal_element_name(s->jibal->elements, s->Z_active));
                return SCRIPT_COMMAND_FAILURE;
            }
            return 0;
        default:
            break;
    }
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Not enough parameters to set detector calibration values.\n", val);
        return SCRIPT_COMMAND_FAILURE;
    }
    c = detector_get_calibration(det, s->Z_active);
    if(!c) {
        jabs_message(MSG_ERROR, "No calibration set for element = %s\n", jibal_element_name(s->jibal->elements, s->Z_active));
        return SCRIPT_COMMAND_FAILURE;
    }
    double value_dbl;
    if(jabs_unit_convert(s->jibal->units, JIBAL_UNIT_TYPE_ANY, argv[0], &value_dbl) < 0) {
        return SCRIPT_COMMAND_FAILURE;
    }
    switch(val) {
        case 's': /* slope */
            if(calibration_set_param(c, CALIBRATION_PARAM_SLOPE, value_dbl)) {
                jabs_message(MSG_ERROR, "Can not set calibration slope.\n");
                return SCRIPT_COMMAND_FAILURE;
            }
            return 1;
        case 'o': /* offset */
            if(calibration_set_param(c, CALIBRATION_PARAM_OFFSET, value_dbl)) {
                jabs_message(MSG_ERROR, "Can not set calibration offset.\n");
                return SCRIPT_COMMAND_FAILURE;
            }
            return 1;
        case 'r': /* resolution */
            if(calibration_set_param(c, CALIBRATION_PARAM_RESOLUTION, value_dbl)) {
                jabs_message(MSG_ERROR, "Can not set calibration resolution.\n");
                return SCRIPT_COMMAND_FAILURE;
            }
            return 1;
        default:
            break;
    }
    jabs_message(MSG_ERROR, "Unhandled value %i in script_set_detector_calibration_val. Report to developer.\n", val);
    return SCRIPT_COMMAND_FAILURE;
}

script_command_status script_set_fit_val(struct script_session *s, int val, int argc, char *const *argv) {
    (void) argv;
    fit_data *fit = s->fit;

    switch(val) {
        case 'n': /* normal */
            fit->phase_start = FIT_PHASE_FAST;
            fit->phase_stop = FIT_PHASE_SLOW;
            return 0;
        case 'f':
            fit->phase_start = FIT_PHASE_FAST;
            fit->phase_stop = FIT_PHASE_FAST;
            return 0;
        case 's':
            fit->phase_start = FIT_PHASE_SLOW;
            fit->phase_stop = FIT_PHASE_SLOW;
            return 0;
        default:
            break;
    }
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Not enough parameters to set fit values.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_set_simulation_val(struct script_session *s, int val, int argc, char *const *argv) {
    (void) argv;
    simulation *sim = s->fit->sim;

    switch(val) {
        case 'd': /* default */
            sim->params = sim_calc_params_defaults(sim->params);
            return 0;
        case 'f': /* fast */
            sim->params = sim_calc_params_defaults_fast(sim->params);
            return 0;
        case 'a': /* accurate */
            sim->params = sim_calc_params_defaults_accurate(sim->params);
            return 0;
        case 'b': /* brisk */
            sim->params = sim_calc_params_defaults_brisk(sim->params);
            return 0;
        case 'i': /* improved */
            sim->params = sim_calc_params_defaults_improved(sim->params);
            return 0;
        default:
            break;
    }
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Not enough parameters to set simulation values.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_set_charge(struct script_session *s, int argc, char *const *argv) {
    (void) argc; /* argc_min set elsewhere */
    double charge;
    if(jabs_unit_convert(s->jibal->units, JIBAL_UNIT_TYPE_CHARGE, argv[0], &charge) < 0) {
        return SCRIPT_COMMAND_FAILURE;
    }
    s->fit->sim->fluence = charge / C_E;
    return 1;
}
script_command_status script_show_var(struct script_session *s, jibal_config_var *var, int argc, char *const *argv) {
    (void) argv;
    (void) argc;
    (void) s;
    DEBUGMSG("Showing var %s of type '%c' (%i)", var->name, var->unit_type, var->unit_type);
    if(var->variable == NULL)
        return SCRIPT_COMMAND_FAILURE;
    switch(var->type) {
        case JIBAL_CONFIG_VAR_NONE:
            break;
        case JIBAL_CONFIG_VAR_PATH:
        case JIBAL_CONFIG_VAR_STRING:
            if(*((void **) var->variable) == NULL)
                break;
            jabs_message(MSG_INFO, "%s = %s\n", var->name, *((char **) var->variable));
            break;
        case JIBAL_CONFIG_VAR_BOOL:
            jabs_message(MSG_INFO, "%s = %s\n", var->name, *((int *) var->variable) ? "true" : "false");
            break;
        case JIBAL_CONFIG_VAR_INT:
            jabs_message(MSG_INFO, "%s = %i\n", var->name, *((int *) var->variable));
            break;
        case JIBAL_CONFIG_VAR_DOUBLE:
            jabs_message(MSG_INFO, "%s = %g\n", var->name, *((double *) var->variable));
            break;
        case JIBAL_CONFIG_VAR_UNIT:
            if(var->unit) {
                jabs_message(MSG_INFO, "%s = %g %s\n", var->name,
                             *((double *) var->variable) / jibal_units_get(s->jibal->units, var->unit_type, var->unit), var->unit);
            } else {
                jabs_message(MSG_INFO, "%s = %g\n", var->name, *((double *) var->variable));
            }
            break;
        case JIBAL_CONFIG_VAR_OPTION:
            jabs_message(MSG_INFO, "%s = %s\n", var->name,
                         jibal_option_get_string(var->option_list, *((int *) var->variable)));
            break;
        case JIBAL_CONFIG_VAR_SIZE:
            jabs_message(MSG_INFO, "%s = %zu\n", var->name, *((size_t *) var->variable));
            break;
    }
    return 0;
}

script_command_status script_enable_var(struct script_session *s, jibal_config_var *var, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    (void) s;
    if(var->type == JIBAL_CONFIG_VAR_BOOL) {
        *((int *) var->variable) = TRUE;
    } else {
        jabs_message(MSG_ERROR, "Variable %s is not boolean.\n", var->name);
    }
    return 0; /* Number of arguments */
}

script_command_status script_disable_var(struct script_session *s, jibal_config_var *var, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    (void) s;
    if(var->type == JIBAL_CONFIG_VAR_BOOL) {
        *((int *) var->variable) = FALSE;
    } else {
        jabs_message(MSG_ERROR, "Variable %s is not boolean.\n", var->name);
    }
    return 0; /* Number of arguments */
}

script_command_status script_execute_command(script_session *s, const char *cmd) {
    int argc = 0;
    script_command_status status;
    DEBUGMSG("Trying to execute command %s", cmd);
    if(!s) {
        jabs_message(MSG_ERROR, "Session has not been initialized.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(jabs_line_is_comment(cmd)) {
        return SCRIPT_COMMAND_SUCCESS;
    }
    char *s_out;
    char **argv = string_to_argv(cmd, &argc, &s_out);
    if(!argv) {
        jabs_message(MSG_ERROR, "Something went wrong in parsing arguments.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc) {
        status = script_execute_command_argv(s, s->commands, argc, argv); /* Note that s->file_depth may be altered (e.g. by script_load_script() */
    } else {
        status = SCRIPT_COMMAND_SUCCESS; /* Doing nothing successfully */
    }
    argv_free(argv, s_out);
    return status;
}

script_command_status script_execute_command_argv(script_session *s, const script_command *commands, int argc, char **argv) {
    if(!s || !commands || !argv)
        return SCRIPT_COMMAND_FAILURE;
    if(argc < 1) {
        jabs_message(MSG_ERROR, "No arguments given.\n");
        return SCRIPT_COMMAND_NOT_FOUND;
    }
    const script_command *cmds = commands;
    const script_command *c_parent = NULL;
    while(argc && cmds) { /* Arguments and subcommands remain. Try to find the right one, if possible. */
        DEBUGMSG("Top level script_execute_command_argv() loop, %i arguments remain (start with %s).", argc, argv[0]);
        if(jabs_line_is_comment(argv[0])) {
            return SCRIPT_COMMAND_SUCCESS;
        }
        const script_command *c = script_command_find(cmds, argv[0]);
        if(!c) {
            DEBUGMSG("Didn't find command %s.", argv[0]);
            script_command_not_found(argv[0], c_parent);
            return SCRIPT_COMMAND_NOT_FOUND;
        }
        while(c) { /* Subcommand found */
            argc--;
            argv++;
            DEBUGMSG("Found command %s.", c->name);
            if(c->f) {
                if(argc < c->argc_min) {
                    jabs_message(MSG_ERROR, "Not enough arguments for command \"%s\", expected minimum of %i\n", c->name, c->argc_min);
                    return EXIT_FAILURE;
                }
                DEBUGMSG("There is a function in command %s. Calling it with %i arguments.", c->name, argc);
                script_command_status status = c->f(s, argc, argv);
                if(status > 0) { /* Positive numbers indicate number of arguments consumed */
                    argc -= status;
                    argv += status;
                }
                DEBUGMSG("Command run, returned %i \"%s\". Number of args remaining: %i",
                         status, script_command_status_to_string(status), argc);
                if(status == SCRIPT_COMMAND_RESET) { /* Special status on script_command_reset(), because now "c" pointer is freed! We must not dereference it anymore. */
                    DEBUGSTR("Return value was SCRIPT_COMMAND_RESET");
                    return SCRIPT_COMMAND_SUCCESS;
                }
                if(status < 0 && status != SCRIPT_COMMAND_NOT_FOUND) { /* Command not found is acceptable, we try to find subcommands later. All other errors cause an immediate return. */
                    DEBUGSTR("Return value is negative, but not SCRIPT_COMMAND_NOT_FOUND");
                    return status;
                }
                if(status == SCRIPT_COMMAND_NOT_FOUND && argc == 0) { /* Command not found, no arguments remain, show an error. */
                    DEBUGSTR("Return value is SCRIPT_COMMAND_NOT_FOUND and no arguments remain (argc == 0)");
                    script_command_not_found(NULL, c);
                    return status;
                }
            } else if(c->var) {
                DEBUGMSG("%s is a var.", c->name);
                if(!c_parent) {
                    jabs_message(MSG_ERROR,
                                 "Command/option \"%s\" is a variable, but there is no parent command at all. This is highly unusual.\n",
                                 c->name);
                    return SCRIPT_COMMAND_FAILURE;
                }
                if(!c_parent->f_var) {
                    jabs_message(MSG_ERROR,
                                 "Command/option \"%s\" is a variable, but there is no function to handle variables in parent command (\"%s\"). This is highly unusual.\n",
                                 c->name, c_parent->name);
                    return SCRIPT_COMMAND_FAILURE;
                }
                DEBUGMSG("Using function in %s. %i arguments. Arguments start with: %s", c_parent->name, argc, argc?argv[0]:"(no arguments)");
                script_command_status status = c_parent->f_var(s, c->var, argc, argv);
                if(status >= 0) { /* Positive numbers indicate number of arguments consumed */
                    DEBUGMSG("Returned positive or zero value, %i arguments consumed.", status);
                    argc -= status;
                    argv += status;
                    DEBUGMSG("Number of arguments remaining: %i.", argc);
                } else {
                    return status;
                }
            } else if(c->val) {
                if(!c_parent) {
                    jabs_message(MSG_ERROR,
                                 "Command/option \"%s\" has a non-zero value, but there is no parent command at all. This is highly unusual.\n",
                                 c->name);
                    return SCRIPT_COMMAND_FAILURE;
                }
                if(!c_parent->f_val) {
                    jabs_message(MSG_ERROR,
                                 "Command/option \"%s\" has a non-zero value, but there is no function to handle values in parent command (\"%s\"). This is highly unusual.\n",
                                 c->name, c_parent->name);
                    return SCRIPT_COMMAND_FAILURE;
                }
                DEBUGMSG("Using function in %s. %i args. Arguments start with: %s", c_parent->name, argc, argc?argv[0]:"(no args)");
                script_command_status status = c_parent->f_val(s, c->val, argc, argv);
                if(status >= 0) { /* Positive numbers indicate number of arguments consumed */
                    argc -= status;
                    argv += status;
                } else {
                    return status;
                }
            } else if(!argc) {
                DEBUGSTR("argc == 0");
                script_command_not_found(NULL, c); /* No function, no nothing, no arguments. */
                return SCRIPT_COMMAND_NOT_FOUND;
            }
            if(c->subcommands) {
                DEBUGSTR("There are subcommands.");
                cmds = c->subcommands;
                c_parent = c;
                break;
            } else {
                DEBUGMSG("there are no subcommands in \"%s\" (not an error). There is val as always: %i.", c->name, c->val);
                c = NULL; /* Moving back to upper level. */
            }
        }
    }
    if(argc) {
        DEBUGMSG("Didn't find command or an error with \"%s\" (%i args remain).", argv[0], argc);
        return SCRIPT_COMMAND_NOT_FOUND;
    }
    return SCRIPT_COMMAND_SUCCESS;
}

const char *script_command_status_to_string(script_command_status status) {
    switch(status) {
        case SCRIPT_COMMAND_NOT_FOUND:
            return "not found";
        case SCRIPT_COMMAND_FAILURE:
            return "failure";
        case SCRIPT_COMMAND_SUCCESS:
            return "success (no arguments consumed)";
        case SCRIPT_COMMAND_EOF:
            return "end-of-file";
        case SCRIPT_COMMAND_EXIT:
            return "exit";
        case SCRIPT_COMMAND_RESET:
            return "reset";
        default:
            break;
    }
    if(status == 1) {
        return "success (1 argument consumed)";
    } else if(status == 2) {
        return "success (2 arguments consumed)";
    } else if(status == 3) {
        return "success (3 arguments consumed)";
    } else if(status > 3) {
        return "success (> 3 arguments consumed)";
    }
    return "unknown";
}

script_command *script_command_new(const char *name, const char *help_text, int val, int argc_min, script_command_status (*f)(struct script_session *, int, char *const *)) {
    if(!name)
        return NULL;
    script_command *c = malloc(sizeof(script_command));
    if(!c)
        return NULL;
    c->name = strdup_non_null(name);
    c->help_text = strdup_non_null(help_text);
    c->f = f;
    c->f_var = NULL;
    c->var = NULL;
    c->val = val;
    c->subcommands = NULL;
    c->next = NULL;
    c->argc_min = argc_min;
    return c;
}

int script_command_set_function(script_command *c, script_command_status (*f)(struct script_session *, int, char *const *)) {
    if(!c)
        return EXIT_FAILURE;
    c->f = f;
    if(c->var) {
        free(c->var);
        c->var = NULL;
    }
    return EXIT_SUCCESS;
}

int script_command_set_var(script_command *c, jibal_config_var_type type, void *variable, const jibal_option *option_list, const char *unit, char unit_type) {
    if(!c->var) {
        c->var = malloc(sizeof(jibal_config_var));
    }
    c->var->variable = variable;
    c->var->type = type;
    c->var->unit = unit;
    c->var->unit_type = unit_type;
    c->var->name = c->name; /* The pointer is shared, so "var" doesn't get its own */
    c->var->option_list = option_list;
    c->f = NULL; /* These guys can't coexist */
    DEBUGVERBOSEMSG("Added command (variable) name %s, type %i, unit %s, unit type '%c' (%i)", c->var->name, c->var->type, c->var->unit, c->var->unit_type, c->var->unit_type);
    return EXIT_SUCCESS;
}

void script_command_free(script_command *c) {
    if(!c)
        return;
    DEBUGVERBOSEMSG("Freeing command \"%s\" (%p)", c->name, (void *)c);
    free(c->var);
    free(c->help_text);
    free(c->name);
    free(c);
}

void script_commands_free(script_command *head) {
    if(!head)
        return;
    struct script_command *stack[SCRIPT_COMMANDS_NESTED_MAX];
    struct script_command *c, *c_old = NULL;
    stack[0] = head;
    size_t i = 0;
    c = stack[0];
    while(c) {
        if(c->subcommands && i < SCRIPT_COMMANDS_NESTED_MAX) { /* Go deeper, push existing pointer to stack */
            stack[i] = c;
            i++;
            c = c->subcommands;
            continue;
        } else {
            c_old = c;
            c = c->next; /* Go to next on the same level */
            script_command_free(c_old);
        }
        while(!c) {
            if(i == 0)
                break;
            i--;
            c = stack[i]->next;
            script_command_free(stack[i]);
        }
    }
}

script_command *script_command_list_find_tail(script_command *head) {
    if(!head)
        return NULL;
    script_command *c = head;
    while(c->next != NULL) {
        c = c->next;
    }
    return c;
}

script_command *script_command_list_merge(script_command *left, script_command *right) {
    script_command *result = NULL;
    script_command *tmp = NULL;

    while(left != NULL && right != NULL) {
        if(strcmp(left->name, right->name) < 0) {
            result = script_command_list_append(result, left);
            tmp = left;
            left = left->next;
            tmp->next = NULL;
        } else {
            result = script_command_list_append(result, right);
            tmp = right;
            right = right->next;
            tmp->next = NULL;
        }
    }
    while(left != NULL) {
        result = script_command_list_append(result, left);
        tmp = left;
        left = left->next;
        tmp->next = NULL;
    }
    while(right != NULL) {
        result = script_command_list_append(result, right);
        tmp = right;
        right = right->next;
        tmp->next = NULL;
    }
    return result;
}

script_command *script_command_list_merge_sort(script_command *head) {
    if(!head)
        return NULL;
    script_command *a[SCRIPT_COMMAND_MERGE_SORT_ARRAY_SIZE];
    memset(a, 0, sizeof(script_command *) * SCRIPT_COMMAND_MERGE_SORT_ARRAY_SIZE);
    script_command *result;
    script_command *next;
    result = head;
    size_t i;
    DEBUGVERBOSEMSG("Sorting command list (head %p = %s)", (void *)head, head->name);
    while(result != NULL) {
        next = result->next;
        result->next = NULL;

        for(i = 0; (i < SCRIPT_COMMAND_MERGE_SORT_ARRAY_SIZE) && (a[i] != NULL); i++) {
            result = script_command_list_merge(a[i], result);
            a[i] = NULL;
        }
        if(i == SCRIPT_COMMAND_MERGE_SORT_ARRAY_SIZE)
            i--;
        a[i] = result;
        result = next;
    }
    result = NULL;
    for(i = 0; i < SCRIPT_COMMAND_MERGE_SORT_ARRAY_SIZE; i++) {
        result = script_command_list_merge(a[i], result);
    }
    return result;
}

script_command *script_command_list_append(script_command *head, script_command *cmd) {
    if(!head) {
        return cmd;
    } else {
        script_command *tail = script_command_list_find_tail(head);
        tail->next = cmd;
    }
    return head;
}

void script_command_list_add_command(script_command **head, script_command *c_new) {
    if(!c_new)
        return;
    if(!c_new->subcommands && !c_new->f && !c_new->var && !c_new->val) { /* Everything needs to be set before calling this function to avoid this warning. */
        DEBUGMSG("Warning: \"%s\" commands/option does not have subcommands and doesn't define a function, variable or a value!", c_new->name);
    }
    if(*head == NULL) {
        *head = c_new;
        return;
    }
    if(c_new->val) { /* Check if vals are multiply defined, because that would suck. */
        for(const script_command *c = *head; c; c = c->next) {
            if(c->val == c_new->val) {
                jabs_message(MSG_WARNING, "Warning: value %i (='%c') defined both in \"%s\" and \"%s\" commands/options!\n", c->val, c->val, c->name, c_new->name);
            }
        }
    }
    script_command *tail = script_command_list_find_tail(*head);
    if(tail) {
        tail->next = c_new;
    }
}

script_command *script_command_list_from_vars_array(const jibal_config_var *vars, jibal_config_var_type type) {
    script_command *head = NULL;
    for(const jibal_config_var *var = vars; var->type != 0; var++) {
        if(type != 0 && var->type != type) { /* Restrict by type */
            continue;
        }
        script_command *c;
        if(var->type == JIBAL_CONFIG_VAR_UNIT) {
            char *help_text;
            if(asprintf(&help_text, "value with unit (type '%c', e.g. %s)", var->unit_type, var->unit) < 0) {
                DEBUGMSG("asprintf failed when adding command %s\n", var->name);
                break;
            }

            c  = script_command_new(var->name, help_text, 0, 0, NULL);
            free(help_text);
        } else if(var->type == JIBAL_CONFIG_VAR_OPTION) {
            if(!var->option_list) {
                DEBUGMSG("Option (%s) without option list!", var->name);
                continue;
            }
            char *help_text;
            if(asprintf(&help_text, "one of:") < 0) {
                DEBUGMSG("asprintf failed when adding command %s\n", var->name);
                break;
            }
            for(const jibal_option *opt = var->option_list; opt->s != NULL; opt++) {
                asprintf_append(&help_text, " \"%s\"", opt->s);
            }
            c  = script_command_new(var->name, help_text, 0, 0, NULL);
            free(help_text);
        } else {
            c = script_command_new(var->name, jibal_config_var_type_name(var->type), 0, 0, NULL);
        }
        script_command_set_var(c, var->type, var->variable, var->option_list, var->unit, var->unit_type);
        script_command_list_add_command(&head, c);
    }
    return head;
}

script_command *script_commands_create(struct script_session *s) {
    fit_data *fit = s->fit;
    simulation *sim = fit->sim;
    script_command *head = NULL;

    script_command *c;
    script_command *c_help = script_command_new("help", "Help.", 0, 0, &script_help);
    script_command_list_add_command(&c_help->subcommands, script_command_new("commands", "List of commands.", 0, 0, &script_help_commands));
    script_command_list_add_command(&c_help->subcommands, script_command_new("version", "Help on (show) version.", 0, 0, &script_help_version));
    script_command_list_add_command(&head, c_help);


#ifdef JABS_PLUGINS
    script_command *c_identify = script_command_new("identify", "Identify something.", 0, 0, NULL);
    script_command_list_add_command(&c_identify->subcommands, script_command_new("plugin", "Identify plugin.", 0, 0, &script_identify_plugin));
    script_command_list_add_command(&head, c_identify);
#endif

    script_command *c_set = script_command_new("set", "Set something.", 0, 0, NULL);
    c_set->f_var = &script_set_var;
    script_command_list_add_command(&c_set->subcommands, script_command_new("aperture", "Set aperture.", 0, 0, &script_set_aperture));
    script_command_list_add_command(&c_set->subcommands, script_command_new("charge", "Set charge.", 0, 1, &script_set_charge)); /* Fluence is fundamental, this is for convenience */
    script_command_list_add_command(&c_set->subcommands, script_command_new("verbosity", "Set verbosity.", 0, 1, &script_set_verbosity));


    script_command *c_channeling = script_command_new("channeling", "Set channeling yield correction (i.e. last layer yield).", 0, 0, NULL);
    script_command_list_add_command(&c_channeling->subcommands, script_command_new("yield", "Set channeling yield (constant).", 0, 0, &script_set_channeling_yield));
    script_command_list_add_command(&c_channeling->subcommands, script_command_new("slope", "Set channeling yield (depth) slope.", 0, 0, &script_set_channeling_slope));
    script_command_list_add_command(&c_set->subcommands, c_channeling);

    script_command *c_detector = script_command_new("detector", "Set detector properties.", 0, 0, &script_set_detector);
    c_detector->f_val = &script_set_detector_val;
    script_command_list_add_command(&c_detector->subcommands, script_command_new("aperture", "Set detector aperture.", 0, 0, &script_set_detector_aperture));
    script_command *c_calibration = script_command_new("calibration", "Set calibration.", 0, 0, &script_set_detector_calibration);
    c_calibration->f_val = &script_set_detector_calibration_val;
    script_command_list_add_command(&c_detector->subcommands, c_calibration);
    script_command_list_add_command(&c_calibration->subcommands, script_command_new(calibration_option[CALIBRATION_LINEAR].s, "Set the calibration to be linear (default).", 'L', 0, NULL));
    script_command_list_add_command(&c_calibration->subcommands, script_command_new("slope", "Set the slope of a linear calibration.", 's', 0, NULL));
    script_command_list_add_command(&c_calibration->subcommands, script_command_new("offset", "Set the offset of a linear calibration.", 'o', 0, NULL));
    script_command_list_add_command(&c_calibration->subcommands, script_command_new("resolution", "Set the resolution.", 'r', 0, NULL));
    script_command_list_add_command(&c_calibration->subcommands, script_command_new(calibration_option[CALIBRATION_POLY].s, "Set the calibration to be a polynomial.", 0, 0, &script_set_detector_calibration_poly));

    script_command_list_add_command(&c_detector->subcommands, script_command_new("column", "Set column number (for data input).", 'c', 0, NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("channels", "Set number of channels.", 'h', 0, NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("compress", "Set compress (summing of channels).", 'C', 0, NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("distance", "Set detector distance from target.", 'd', 0, NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("foil", "Set detector foil.", 0, 0, &script_set_detector_foil));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("type", "Set detector type.", 't', 0, NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("slope", "Set detector calibration slope.", 'S', 0, NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("offset", "Set detector calibration offset.", 'O', 0, NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("resolution", "Set detector resolution (FWHM).", 'r', 0, NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("solid", "Set detector solid angle.", 's', 0, NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("theta", "Set detector (scattering) angle.", 'T', 0, NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("length", "Set detector length (for ToF).", 'l', 0, NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("beta", "Set exit angle (angle of ion in sample) for this detector (manually).", 'b', 0, NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("phi", "Set detector azimuth angle, 0 = IBM, 90 deg = Cornell.", 'p', 0, NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("name", "Set detector name.", 0, 0, &script_set_detector_name));
    script_command_list_add_command(&c_set->subcommands, c_detector);

    script_command_list_add_command(&c_set->subcommands, script_command_new("ion", "Set incident ion (isotope).", 0, 0, &script_set_ion));

    script_command *c_set_fit = script_command_new("fit", "Set fit related things.", 0, 0, NULL);
    c_set_fit->f_val = &script_set_fit_val;
    script_command_list_add_command(&c_set_fit->subcommands, script_command_new("normal", "Normal two-phase fitting.", 'n', 0, NULL));
    script_command_list_add_command(&c_set_fit->subcommands, script_command_new("slow", "One phase fitting (slow phase only).", 's', 0, NULL));
    script_command_list_add_command(&c_set_fit->subcommands, script_command_new("fast", "One phase fitting (fast phase only).", 'f', 0, NULL));
    script_command_list_add_command(&c_set->subcommands, c_set_fit);

    script_command_list_add_command(&c_set->subcommands, script_command_new("sample", "Set sample.", 0, 0, &script_set_sample));

    script_command *c_set_simulation = script_command_new("simulation", "Set simulation related things.", 0, 0, NULL);
    c_set_simulation->f_val = &script_set_simulation_val;
    script_command_list_add_command(&c_set_simulation->subcommands, script_command_new("defaults", "Set default calculation parameters.", 'd', 0, NULL));
    script_command_list_add_command(&c_set_simulation->subcommands, script_command_new("brisk", "Set slightly faster calculation parameters.", 'b', 0, NULL));
    script_command_list_add_command(&c_set_simulation->subcommands, script_command_new("fast", "Set fast calculation parameters.", 'f', 0, NULL));
    script_command_list_add_command(&c_set_simulation->subcommands, script_command_new("accurate", "Set the most accurate calculation parameters.", 'a', 0, NULL));
    script_command_list_add_command(&c_set_simulation->subcommands, script_command_new("improved", "Set more accurate calculation parameters.", 'i', 0, NULL));
    script_command_list_add_command(&c_set->subcommands, c_set_simulation);

    script_command_list_add_command(&c_set->subcommands, script_command_new("stopping", "Set (assign) stopping or straggling.", 0, 0, &script_set_stopping));

    const jibal_config_var vars[] = {
            {JIBAL_CONFIG_VAR_DOUBLE, "fluence",                       0,     0,                               &sim->fluence,                               NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "energy",                        "keV", JIBAL_UNIT_TYPE_ENERGY,          &sim->beam_E,                                NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "energy_broad",                  "keV", JIBAL_UNIT_TYPE_ENERGY,          &sim->beam_E_broad,                          NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "emin",                          "keV", JIBAL_UNIT_TYPE_ENERGY,          &sim->emin,                                  NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "alpha",                         "deg", JIBAL_UNIT_TYPE_ANGLE,           &sim->sample_theta,                          NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "phi",                           "deg", JIBAL_UNIT_TYPE_ANGLE,           &sim->sample_phi,                            NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "erd",                           0,     0,                               &sim->erd,                                   NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "rbs",                           0,     0,                               &sim->rbs,                                   NULL},
            {JIBAL_CONFIG_VAR_SIZE,   "maxiter",                       0,     0,                               &fit->n_iters_max,                           NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "xtolerance",                    0,     0,                               &fit->xtol,                                  NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "chisq_tolerance",               0,     0,                               &fit->chisq_tol,                             NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "chisq_fast_tolerance",          0,     0,                               &fit->chisq_fast_tol,                        NULL},
            {JIBAL_CONFIG_VAR_SIZE,   "n_bricks_max",                  0,     0,                               &sim->params->n_bricks_max,                  NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "ds",                            0,     0,                               &sim->params->ds,                            NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "rk4",                           0,     0,                               &sim->params->rk4,                           NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "brick_width_sigmas",            0,     0,                               &sim->params->brick_width_sigmas,            NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "sigmas_cutoff",                 0,     0,                               &sim->params->sigmas_cutoff,                 NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "incident_stop_step",            "keV", JIBAL_UNIT_TYPE_ENERGY,          &sim->params->incident_stop_params.step,     NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "incident_stop_step_sigmas",     0,     0,                               &sim->params->incident_stop_params.sigmas,   NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "incident_stop_step_min",        "keV", JIBAL_UNIT_TYPE_ENERGY,          &sim->params->incident_stop_params.min,      NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "incident_stop_step_max",        "keV", JIBAL_UNIT_TYPE_ENERGY,          &sim->params->incident_stop_params.max,      NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "exiting_stop_step",             "keV", JIBAL_UNIT_TYPE_ENERGY,          &sim->params->exiting_stop_params.step,      NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "exiting_stop_step_sigmas",      0,     0,                               &sim->params->exiting_stop_params.sigmas,    NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "exiting_stop_step_min",         "keV", JIBAL_UNIT_TYPE_ENERGY,          &sim->params->exiting_stop_params.min,       NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "exiting_stop_step_max",         "keV", JIBAL_UNIT_TYPE_ENERGY,          &sim->params->exiting_stop_params.max,       NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "ds_incident_stop_step_factor",  0,     0,                               &sim->params->ds_incident_stop_step_factor,  NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "nuclear_stopping_accurate",     0,     0,                               &sim->params->nuclear_stopping_accurate,     NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "mean_conc_and_energy",          0,     0,                               &sim->params->mean_conc_and_energy,          NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "geostragg",                     0,     0,                               &sim->params->geostragg,                     NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "beta_manual",                   "deg", JIBAL_UNIT_TYPE_ANGLE,           &sim->params->beta_manual,                   NULL},
            {JIBAL_CONFIG_VAR_SIZE,   "cs_n_stragg_steps",             0,     0,                               &sim->params->cs_n_stragg_steps,             NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "gaussian_accurate",             0,     0,                               &sim->params->gaussian_accurate,             NULL},
            {JIBAL_CONFIG_VAR_SIZE,   "int_cs_max_intervals",          0,     0,                               &sim->params->int_cs_max_intervals,          NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "int_cs_accuracy",               0,     0,                               &sim->params->int_cs_accuracy,               NULL},
            {JIBAL_CONFIG_VAR_SIZE,   "int_cs_stragg_max_intervals",   0,     0,                               &sim->params->int_cs_stragg_max_intervals,   NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "int_cs_stragg_accuracy",        0,     0,                               &sim->params->int_cs_stragg_accuracy,        NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "cs_adaptive",                   0,     0,                               &sim->params->cs_adaptive,                   NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "cs_energy_step_max",            "keV", JIBAL_UNIT_TYPE_ENERGY,          &sim->params->cs_energy_step_max,            NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "cs_depth_step_max",             "tfu", JIBAL_UNIT_TYPE_LAYER_THICKNESS, &sim->params->cs_depth_step_max,             NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "cs_stragg_step_sigmas",         0,     0,                               &sim->params->cs_stragg_step_sigmas,         NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "reaction_file_angle_tolerance", "deg", JIBAL_UNIT_TYPE_ANGLE,           &sim->params->reaction_file_angle_tolerance, NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "bricks_skip_zero_conc_ranges",  0,     0,                               &sim->params->bricks_skip_zero_conc_ranges,  NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "screening_tables",              0,     0,                               &sim->params->screening_tables,              NULL},
            {JIBAL_CONFIG_VAR_OPTION, "cs_rbs",                        0,     0,                               &sim->cs_rbs,                                jabs_cs_types},
            {JIBAL_CONFIG_VAR_OPTION, "cs_erd",                        0,     0,                               &sim->cs_erd,                                jabs_cs_types},
            {JIBAL_CONFIG_VAR_NONE, NULL,                              0,     0, NULL,                                                                      NULL}
    };
    c = script_command_list_from_vars_array(vars, 0);
    script_command_list_add_command(&c_set->subcommands, c);
    script_command_list_add_command(&head, c_set); /* End of "set" commands */

    script_command *c_load = script_command_new("load", "Load something.", 0, 0, NULL);
    script_command_list_add_command(&c_load->subcommands, script_command_new("experimental", "Load an experimental spectrum.", 0, 0, &script_load_experimental));
    script_command_list_add_command(&c_load->subcommands, script_command_new("reference", "Load a reference spectrum.", 0, 0, &script_load_reference));
    script_command_list_add_command(&c_load->subcommands, script_command_new("script", "Load (run) a script.", 0, 0, &script_load_script));
    script_command_list_add_command(&c_load->subcommands, script_command_new("sample", "Load a sample.", 0, 0, &script_load_sample));
    script_command *c_reaction = script_command_new("reaction", "Load a reaction from R33 file.", 0, 0, &script_load_reaction);
    script_command_list_add_command(&c_load->subcommands, c_reaction);
#ifdef JABS_PLUGINS
    script_command_list_add_command(&c_reaction->subcommands, script_command_new("plugin", "Load a reaction from a plugin.", 0, 0, &script_load_reaction_plugin));
#endif
    script_command_list_add_command(&c_load->subcommands, script_command_new("roughness", "Load layer thickness table (roughness) from a file.", 0, 0, &script_load_roughness));
    script_command_list_add_command(&head, c_load); /* End of "load" commands */

    script_command *c_show = script_command_new("show", "Show information on things.", 0, 0, NULL);
    script_command_list_add_command(&c_show->subcommands, script_command_new("aperture", "Show aperture.", 0, 0, &script_show_aperture));
    script_command_list_add_command(&c_show->subcommands, script_command_new("calc_params", "Show calculation parameters.", 0, 0, &script_show_calc_params));
    script_command_list_add_command(&c_show->subcommands, script_command_new("detector", "Show detector.", 0, 0, &script_show_detector));

    script_command *c_fit = script_command_new("fit", "Show fit results.", 0, 0, &script_show_fit);
    script_command_list_add_command(&c_show->subcommands, c_fit);
    script_command_list_add_command(&c_fit->subcommands, script_command_new("variables", "Show possible fit variables.", 0, 0, &script_show_fit_variables));
    script_command_list_add_command(&c_fit->subcommands, script_command_new("ranges", "Show fit ranges.", 0, 0, &script_show_fit_ranges));

    script_command_list_add_command(&c_show->subcommands, script_command_new("reactions", "Show reactions.", 0, 0, &script_show_reactions));

    script_command *c_show_sample = script_command_new("sample", "Show sample.", 0, 0, &script_show_sample);
    script_command_list_add_command(&c_show_sample->subcommands, script_command_new("profile", "Show sample profile.", 0, 0, &script_show_sample_profile));
    script_command_list_add_command(&c_show->subcommands, c_show_sample);

    script_command_list_add_command(&c_show->subcommands, script_command_new("simulation", "Show simulation.", 0, 0, &script_show_simulation));
    script_command_list_add_command(&c_show->subcommands, script_command_new("stopping", "Show stopping (GSTO) assignments.", 0, 0, &script_show_stopping));

    script_command *c_show_variable = script_command_new("variable", "Show variable.", 0, 0, NULL);
    c_show_variable->f_var = script_show_var;
    c = script_command_list_from_vars_array(vars, 0);
    script_command_list_add_command(&c_show_variable->subcommands, c);
    script_command_list_add_command(&c_show->subcommands, c_show_variable);
    script_command_list_add_command(&head, c_show); /* End of "show" commands */

    script_command_list_add_command(&head, script_command_new("exit", "Exit.", 0, 0, &script_exit));

    c = script_command_new("save", "Save something.", 0, 0, NULL);
    script_command_list_add_command(&c->subcommands, script_command_new("calibrations", "Save detector calibrations.", 0, 0, &script_save_calibrations));
    script_command_list_add_command(&c->subcommands, script_command_new("sample", "Save sample.", 0, 0, &script_save_sample));
    script_command_list_add_command(&c->subcommands, script_command_new("simulation", "Save simulation.", 0, 0, &script_save_simulation));
    script_command_list_add_command(&c->subcommands, script_command_new("spectra", "Save spectra.", 0, 0, &script_save_spectra));
    script_command_list_add_command(&head, c); /* End of "save" commands */

    c = script_command_new("test", "Test something.", 0, 0, NULL);
    script_command_list_add_command(&c->subcommands, script_command_new("reference", "Test simulated spectrum against reference spectrum.", 0, 0, &script_test_reference));
    script_command_list_add_command(&c->subcommands, script_command_new("roi", "Test ROI.", 0, 0, &script_test_roi));
    script_command_list_add_command(&head, c); /* End of "test" commands */

    c = script_command_new("add", "Add something.", 0, 0, NULL);
    script_command *c_add_detector = script_command_new("detector", "Add a detector.", 0, 0, NULL);
    script_command_list_add_command(&c_add_detector->subcommands, script_command_new("default", "Add a default detector.", 0, 0, &script_add_detector_default));
    script_command_list_add_command(&c->subcommands, c_add_detector);


    script_command *c_add_fit = script_command_new("fit", "Add something related to fit.", 0, 0, NULL);
    script_command_list_add_command(&c_add_fit->subcommands, script_command_new("range", "Add a fit range", 0, 0, &script_add_fit_range));
    script_command_list_add_command(&c->subcommands, c_add_fit);

    script_command_list_add_command(&c->subcommands, script_command_new("reaction", "Add a reaction.", 0, 0, &script_add_reaction));
    script_command_list_add_command(&c->subcommands, script_command_new("reactions", "Add reactions (of some type).", 0, 0, &script_add_reactions));
    script_command_list_add_command(&head, c); /* End of "add" commands */

    c = script_command_new("remove", "Remove something.", 0, 0, NULL);
    script_command_list_add_command(&c->subcommands, script_command_new("reaction", "Remove reaction.", 0, 0, &script_remove_reaction));
    script_command_list_add_command(&head, c);

    c = script_command_new("reset", "Reset something (or everything).", 0, 0, &script_reset);
    script_command_list_add_command(&c->subcommands, script_command_new("detectors", "Reset detectors.", 0, 0, &script_reset_detectors));
    script_command_list_add_command(&c->subcommands, script_command_new("experimental", "Reset experimental spectra.", 0, 0, &script_reset_experimental));
    script_command_list_add_command(&c->subcommands, script_command_new("reference", "Reset reference spectrum.", 0, 0, &script_reset_reference));
    script_command_list_add_command(&c->subcommands, script_command_new("fit", "Reset fit (ranges).", 0, 0, &script_reset_fit));
    script_command_list_add_command(&c->subcommands, script_command_new("reactions", "Reset reactions.", 0, 0, &script_reset_reactions));
    script_command_list_add_command(&c->subcommands, script_command_new("sample", "Reset sample.", 0, 0, &script_reset_sample));
    script_command_list_add_command(&c->subcommands, script_command_new("stopping", "Reset stopping assignments.", 0, 0, &script_reset_stopping));
    script_command_list_add_command(&head, c);

    c = script_command_new("fit", "Do a fit.", 0, 0, script_fit);
    script_command_list_add_command(&head, c);

    script_command *c_enable = script_command_new("enable", "Set boolean variable to true.", 0, 0, NULL);
    c_enable->f_var = &script_enable_var;
    c = script_command_list_from_vars_array(vars, JIBAL_CONFIG_VAR_BOOL);
    script_command_list_add_command(&c_enable->subcommands, c);
    script_command_list_add_command(&head, c_enable);

    script_command *c_disable = script_command_new("disable", "Set boolean variable to true.", 0, 0, NULL);
    c_disable->f_var = &script_disable_var;
    c = script_command_list_from_vars_array(vars, JIBAL_CONFIG_VAR_BOOL);
    script_command_list_add_command(&c_disable->subcommands, c);
    script_command_list_add_command(&head, c_disable);

    c = script_command_new("roi", "Show information from a region of interest.", 0, 0, script_roi);
    script_command_list_add_command(&head, c);

    c = script_command_new("simulate", "Run a simulation.", 0, 0, script_simulate);
    script_command_list_add_command(&head, c);

    script_command *c_split = script_command_new("split", "Split something.", 0, 0, NULL);
    script_command *c_split_sample = script_command_new("sample", "Split something sample related.", 0, 0, NULL);
    script_command_list_add_command(&c_split_sample->subcommands, script_command_new("elements", "Split materials down to their constituent elements.", 0, 0, script_split_sample_elements));
    script_command_list_add_command(&c_split->subcommands, c_split_sample);
    script_command_list_add_command(&head, c_split);

    script_command_list_add_command(&head, script_command_new("cwd", "Display current working directory.", 0, 0, script_cwd));
    script_command_list_add_command(&head, script_command_new("pwd", "Display current working directory.", 0, 0, script_cwd));
    script_command_list_add_command(&head, script_command_new("cd", "Change current working directory.", 0, 0, script_cd));
    script_command_list_add_command(&head, script_command_new("idf2jbs", "Convert IDF file to a JaBS script.", 0, 0, script_idf2jbs));
    script_command_list_add_command(&head, script_command_new("kinematics", "Calculate kinematics.", 0, 0, script_kinematics));
    return script_commands_sort_all(head);
}

script_command *script_commands_sort_all(script_command *head) {
    if(!head)
        return NULL;
    struct script_command *stack[SCRIPT_COMMANDS_NESTED_MAX];
    struct script_command *c;
    head = script_command_list_merge_sort(head);
    stack[0] = head;
    size_t i = 0;
    c = stack[0];
    while(c) {
        if(c->subcommands && i < SCRIPT_COMMANDS_NESTED_MAX) { /* Go deeper, push existing pointer to stack */
            stack[i] = c;
            i++;
            DEBUGVERBOSEMSG("Sorting all commands, level %zu, subcommands of %s", i, c->name);
            c->subcommands = script_command_list_merge_sort(c->subcommands);
            c = c->subcommands;
            continue;
        } else {
            c = c->next; /* Go to next on the same level */
        }
        while(!c) {
            if(i == 0)
                break;
            i--;
            c = stack[i]->next;
        }
    }
    return head;
}

void script_commands_print(const struct script_command *commands) {
    if(!commands)
        return;
    size_t len_max = 0;
    for(const struct script_command *c = commands; c; c = c->next) {
        if(!c->help_text)
            continue;
        len_max = JABS_MAX(len_max, strlen(c->name));
    }
    for(const struct script_command *c = commands; c; c = c->next) {
        if(!c->help_text)
            continue;
        jabs_message(MSG_INFO, " %*s    %s\n", len_max, c->name, c->help_text);
    }
}

void script_print_command_tree(const struct script_command *commands) {
    const struct script_command *stack[SCRIPT_COMMANDS_NESTED_MAX];
    const struct script_command *c;
    stack[0] = commands;
    size_t i = 0;
    c = stack[0];
    while(c) {
        if(c->f || c->var || c->val) { /* If none of these is set, we shouldn't print the command name at all */
            for(size_t j = 0; j < i; j++) {
                jabs_message(MSG_INFO, "%s ", stack[j]->name);
            }
            jabs_message(MSG_INFO, "%s\n", c->name);
        } else {
#ifdef DEBUG
            DEBUGSTR("The heck is this?");
            abort();
#endif
        }
        if(c->subcommands && i < SCRIPT_COMMANDS_NESTED_MAX) { /* Go deeper, push existing pointer to stack */
            stack[i] = c;
            i++;
            c = c->subcommands;
            continue;
        } else {
            c = c->next; /* Go to next on the same level */
        }
        while(!c) {
            if(i == 0)
                break;
            i--;
            c = stack[i]->next;
        }
    }
}

script_command_status script_load_script(script_session *s, int argc, char *const *argv) {
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Usage: load script [file]\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(script_session_load_script(s, argv[0])) {
        return SCRIPT_COMMAND_FAILURE;
    }
    return 1;
}

script_command_status script_load_sample(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Usage: load sample [file]\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    sample_model *sm = sample_model_from_file(fit->jibal, argv[0]);
    if(!sm) {
        jabs_message(MSG_INFO, "Sample load from \"%s\" failed.\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    if(fit_data_set_sample_model(fit, sm)) {
        jabs_message(MSG_WARNING, "Setting sample fails due to sample model sanity check.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    return 1;
}

script_command_status script_load_experimental(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    const int argc_orig = argc;
    size_t i_det = 0;
    if(script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc < 1) {
        jabs_message(MSG_ERROR,
                     "Usage: load experimental {<detector>} <filename>\nIf no detector number is given, experimental data is loaded for all detectors from the same file.\nUse 'set detector column' to control which columns are read.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc != argc_orig) { /* Detector number given */
        if(fit_data_load_exp(s->fit, i_det, argv[0])) {
            return SCRIPT_COMMAND_FAILURE;
        }
    } else { /* Load all detectors from same file (with multiple columns). TODO: make a better multicolumn data reader that can do this in one pass. */
        for(i_det = 0; i_det < fit->sim->n_det; i_det++) {
            if(fit_data_load_exp(s->fit, i_det, argv[0])) {
                return SCRIPT_COMMAND_FAILURE;
            } else {
                const jabs_histogram *h = fit->exp[i_det];
                jabs_message(MSG_VERBOSE, "Detector %zu: experimental spectrum with %zu channels loaded.\n", i_det + 1, h?h->n:0);
                calibration_apply_to_histogram(detector_get_calibration(sim_det(fit->sim, i_det), JIBAL_ANY_Z), fit_data_exp(fit, i_det));
            }
        }
    }
    argc--;
    argv++;
    return argc_orig - argc; /* Number of arguments */
}

script_command_status script_load_reference(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    const int argc_orig = argc;
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Usage: load reference <filename>\n");
    }
    const char *filename = argv[0];
    jabs_histogram *h =spectrum_read(filename, 0, CHANNELS_MAX_DEFAULT, 1, 1);
    if(!h) {
        jabs_message(MSG_ERROR, "Reading reference spectrum from file \"%s\" was not successful.\n", filename);
        return EXIT_FAILURE;
    }
    jabs_message(MSG_VERBOSE, "Reference spectrum from file \"%s\" loaded. %zu channels.\n", filename, h->n);
    jabs_histogram_free(fit->ref);
    fit->ref = h;
    argc--;
    argv++;
    return argc_orig - argc; /* Number of arguments */
}

script_command_status script_load_reaction(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    const int argc_orig = argc;
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Usage: load reaction <file>\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(strcmp(argv[0], "plugin") == 0) {
        return 0;
    }
    const char *filename = argv[0];
    argc--;
    argv++;
    r33_file *rfile = r33_file_read(filename);
    if(!rfile) {
        jabs_message(MSG_ERROR, "Could not load R33 from file \"%s\".\n", filename);
        return SCRIPT_COMMAND_FAILURE;
    }
    reaction *r = r33_file_to_reaction(s->jibal->isotopes, rfile);
    r33_file_free(rfile);
    if(!r) {
        jabs_message(MSG_ERROR, "Could not convert R33 file to a reaction!\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    jabs_message(MSG_INFO, "File: %s has a reaction with %s -> %s, product %s, theta %g deg\n",
                 filename, r->incident->name, r->target->name, r->product->name, r->theta / C_DEG);
    if(reaction_modifiers_from_argv(s->jibal, r, &argc, &argv)) {
        reaction_free(r);
        jabs_message(MSG_ERROR, "Could not parse reaction modifiers.");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(sim_reactions_add_reaction(fit->sim, r, FALSE)) {
        reaction_free(r);
        jabs_message(MSG_ERROR, "Could not parse reaction modifiers.");
        return SCRIPT_COMMAND_FAILURE;
    }
    return argc_orig - argc; /* Number of arguments */
}

script_command_status script_load_roughness(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    sample_model *sm = fit->sm;
    if(argc < 2) {
        jabs_message(MSG_ERROR, "Usage: load roughness <layer> <file>\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(!sm) {
        jabs_message(MSG_ERROR, "No sample has been set, can't use this command yet.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    size_t n = strtoull(argv[0], NULL, 10);
    if(n == 0 || n > sm->n_ranges) {
        jabs_message(MSG_ERROR, "Layer number must be between 1 and %zu (number of ranges in current sample).\n", sm->n_ranges);
        return SCRIPT_COMMAND_FAILURE;
    }
    sample_range *range = &(sm->ranges[n - 1]);
    if(roughness_set_from_file(&range->rough, argv[1])) {
        jabs_message(MSG_ERROR, "Setting roughness from file failed.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    range->x = thickness_probability_table_areal_density(range->rough.file->tpd); /* Update thickness of layer to correspond the average areal density */
    if(range->rough.file && range->rough.file->tpd) {
        jabs_message(MSG_INFO, "Layer %zu: roughness from file \"%s\", containing %zu data points, loaded.\n", n, range->rough.file->filename, range->rough.file->tpd->n);
        const thick_prob_dist *tpd = range->rough.file->tpd;
        jabs_message(MSG_INFO, "  i | thickness (tfu) | weight(%%)\n");
        for(size_t i = 0; i < tpd->n; i++) {
            jabs_message(MSG_INFO, "%3zu | %15.3lf | %9.3lf\n", i + 1, tpd->p[i].x / C_TFU, tpd->p[i].prob * 100.0);
        }
        jabs_message(MSG_INFO, "Average areal density: %g tfu\n", range->x / C_TFU);
    }
    return 2;
}

script_command_status script_reset_reactions(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    sim_reactions_free(s->fit->sim);
    return 0;
}

script_command_status script_reset_detectors(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    for(size_t i_det = 0; i_det < s->fit->sim->n_det; i_det++) {
        detector_free(sim_det(s->fit->sim, i_det));
    }
    s->fit->sim->n_det = 0;
    return 0;
}

script_command_status script_reset_fit(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    fit_data_reset(s->fit);
    return 0;
}

script_command_status script_reset_sample(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    struct fit_data *fit = s->fit;
    sample_model_free(fit->sm);
    fit->sm = NULL;
    return 0;
}

script_command_status script_reset_stopping(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    jibal_gsto_assign_clear_all(s->fit->jibal->gsto);
    return 0;
}

script_command_status script_reset_experimental(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    fit_data_exp_reset(s->fit);
    return 0;
}

script_command_status script_reset_reference(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    jabs_histogram_free(s->fit->ref);
    s->fit->ref = NULL;
    return 0;
}

script_command_status script_reset(script_session *s, int argc, char *const *argv) {
    if(!s) {
        return SCRIPT_COMMAND_FAILURE;
    }
    struct fit_data *fit = s->fit;
    (void) argc; /* Unused */
    (void) argv; /* Unused */
    if(argc > 0) {
        return SCRIPT_COMMAND_NOT_FOUND;
    }
    if(!fit) {
        return SCRIPT_COMMAND_FAILURE;
    }
    DEBUGSTR("Resetting everything!\n");
    fit_data_exp_reset(fit);
    jabs_histogram_free(fit->ref);
    fit->ref = NULL;

    sample_model_free(fit->sm);
    fit->sm = NULL;
    sim_free(fit->sim);
    fit_data_free(s->fit);
    s->fit = fit_data_new(s->jibal, sim_init(s->jibal));
    if(!s->fit) {
        jabs_message(MSG_ERROR, "Reset fails due to an allocation issue. This should not happen.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    script_commands_free(s->commands);
    s->commands = script_commands_create(s);
    jibal_gsto_assign_clear_all(s->jibal->gsto);
    return SCRIPT_COMMAND_RESET;
}

script_command_status script_show_sample(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    if(argc > 0) {
        return SCRIPT_COMMAND_NOT_FOUND;
    }
    struct fit_data *fit = s->fit;
    if(!fit->sm) {
        jabs_message(MSG_WARNING, "No sample has been set.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
#ifdef DEBUG
    char *sample_str = sample_model_to_string(fit->sm);
        fprintf(stderr, "Sample: %s\n", sample_str);
        free(sample_str);
#endif
    jabs_message(MSG_INFO, "Sample model (use \"save sample\" to save):\n");
    if(sample_model_print(NULL, fit->sm, MSG_INFO)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    if(fit->sm->type == SAMPLE_MODEL_LAYERED) {
        jabs_message(MSG_INFO, "\nLayers:\n");
        sample *sample = sample_from_sample_model(fit->sm);
        sample_print_thicknesses(NULL, sample, MSG_INFO);
        sample_free(sample);
    }
    return 0;
}

script_command_status script_show_sample_profile(struct script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    struct fit_data *fit = s->fit;
    if(!fit->sm) {
        jabs_message(MSG_WARNING, "No sample has been set.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    sample *sample = sample_from_sample_model(fit->sm);
    sample_print(sample, FALSE, MSG_INFO);
    sample_free(sample);
    return 0;
}

script_command_status script_show_simulation(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    sim_print(s->fit->sim, MSG_INFO);
    return 0;
}

script_command_status script_show_stopping(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    jibal_gsto *gsto = s->jibal->gsto;
    int n = 0;
    jabs_message(MSG_INFO, "List of assigned stopping and straggling files:\n");
    for(int Z1 = 1; Z1 <= gsto->Z1_max; Z1++) {
        for(int Z2 = 1; Z2 <= gsto->Z2_max; Z2++) {
            gsto_file_t *file_sto = jibal_gsto_get_assigned_file(gsto, GSTO_STO_ELE, Z1, Z2);
            gsto_file_t *file_stg = jibal_gsto_get_assigned_file(gsto, GSTO_STO_STRAGG, Z1, Z2);
            if(!file_sto && !file_stg) {
                continue; /* Nothing to do, nothing assigned */
            }
            jabs_message(MSG_INFO, "  Z1=%i (%s), Z2=%i (%s): ", Z1, jibal_element_name(gsto->elements, Z1), Z2,
                         jibal_element_name(gsto->elements, Z2));
            if(file_sto) {
                jabs_message(MSG_INFO, "Stopping file %s.", file_sto->name);
                n++;
            }
            if(file_stg) {
                jabs_message(MSG_INFO, "%sStraggling file %s.", file_sto ? " " : "", file_stg->name);
                n++;
            }
            jabs_message(MSG_INFO, "\n");
        }
    }
    jabs_message(MSG_INFO, "Total of %i assignments.\n", n);
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_show_fit(script_session *s, int argc, char *const *argv) {
    (void) argv;
    if(argc == 0) {
        if(!s->fit->fit_params || s->fit->fit_params->n_active == 0) {
            jabs_message(MSG_INFO, "No fit results to show.\n");
        } else {
            fit_params_print_final(s->fit->fit_params);
        }
    }
    return 0;
}

script_command_status script_show_fit_variables(script_session *s, int argc, char *const *argv) {
    const int argc_orig = argc;
    fit_params *p_all = fit_params_all(s->fit);
    char *pattern = NULL;
    if(argc > 0) {
        pattern = argv[0];
        argv++;
        argc--;
    }
    fit_params_print(p_all, FALSE, pattern, MSG_INFO);
    fit_params_free(p_all);
    return argc_orig - argc;
}

script_command_status script_show_fit_ranges(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    fit_data_print(s->fit, MSG_INFO);
    return 0;
}

script_command_status script_show_aperture(struct script_session *s, int argc, char *const *argv) {
    (void) argv;
    struct fit_data *fit = s->fit;
    const int argc_orig = argc;
    char *aperture_str = aperture_to_string(fit->sim->beam_aperture);
    jabs_message(MSG_INFO, "aperture %s\n", aperture_str);
    free(aperture_str);
    return argc_orig - argc;
}

script_command_status script_show_calc_params(script_session *s, int argc, char * const *argv) {
    (void) argv;
    const int argc_orig = argc;
    sim_calc_params_print(s->fit->sim->params, MSG_INFO);
    return argc_orig - argc;
}


script_command_status script_show_detector(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    size_t i_det = 0;
    const int argc_orig = argc;
    if(script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    if(fit->sim->n_det == 0) {
        jabs_message(MSG_INFO, "No detectors have been defined.\n");
    } else if(argc != argc_orig || fit->sim->n_det == 1) { /* Show information on particular detector if a number was given or if we only have one detector. */
        if(detector_print(s->jibal, sim_det(fit->sim, i_det))) {
            jabs_message(MSG_ERROR, "No detectors set or other error.\n");
        }
    } else {
        jabs_message(MSG_INFO, "  # |         name | col | theta |   phi |   solid |     type | resolution | calibration\n");
        jabs_message(MSG_INFO, "    |              |     |   deg |   deg |     msr |          |            |            \n");
        for(size_t i = 0; i < fit->sim->n_det; i++) {
            detector *det = sim_det(fit->sim, i);
            if(!det)
                continue;
            char *calib_str = calibration_to_string(det->calibration);
            char *reso_str = detector_resolution_to_string(det, JIBAL_ANY_Z);
            jabs_message(MSG_INFO, "%3zu | %12s | %3zu | %5.1lf | %5.1lf | %7.3lf | %8s | %10s | %s\n",
                         i + 1, det->name, det->column, det->theta / C_DEG, det->phi / C_DEG, det->solid / C_MSR, detector_type_name(det), reso_str, calib_str);
            free(calib_str);
            free(reso_str);
        }
        jabs_message(MSG_INFO, "Use 'show detector <number>' to get more information on a particular detector or 'reset detectors' to remove all of them.\n");
    }
    return argc_orig - argc; /* Number of arguments */
}

script_command_status script_show_reactions(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    struct fit_data *fit = s->fit;
    if(fit->sim->n_reactions == 0) {
        jabs_message(MSG_INFO, "No reactions.\n");
    }
    reactions_print(fit->sim->reactions, fit->sim->n_reactions);
    return 0;
}

script_command_status script_set_ion(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Usage: set ion <isotope>\nExample: set ion 4He\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    const jibal_isotope *isotope = jibal_isotope_find(fit->jibal->isotopes, argv[0], 0, 0);
    if(!isotope) {
        jabs_message(MSG_ERROR, "No such isotope: %s\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }

    if(s->fit->sim->beam_isotope != isotope) {
        s->fit->sim->beam_isotope = isotope;
        if(s->fit->sim->n_reactions > 0) {
            jabs_message(MSG_WARNING, "Reactions were reset automatically, since the ion was changed.\n");
            sim_reactions_free(fit->sim);
        }
    }

    return 1;
}

script_command_status script_set_channeling_yield(script_session *s, int argc, char *const *argv) {
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Usage: set channeling yield <yield>\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    sample_model *sm = s->fit->sm;
    if(!sm || sm->n_ranges == 0) {
        jabs_message(MSG_ERROR, "Sample needs to be set first.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    sm->ranges[sm->n_ranges - 1].yield = strtod(argv[0], NULL);
    return 1;
}

script_command_status script_set_channeling_slope(script_session *s, int argc, char *const *argv) {
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Usage: set channeling slope <slope>\nSlope is assumed to have units of 1/tfu, no unit can be provided.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    sample_model *sm = s->fit->sm;
    if(!sm || sm->n_ranges == 0) {
        jabs_message(MSG_ERROR, "Sample needs to be set first.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    sm->ranges[sm->n_ranges - 1].yield_slope = strtod(argv[0], NULL);
    return 1;
}

script_command_status script_set_aperture(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    const int argc_orig = argc;
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Usage: set aperture <type> {width|height|diameter} ...\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    fit->sim->beam_aperture = aperture_set_from_argv(s->jibal, fit->sim->beam_aperture, &argc, &argv);
    if(!fit->sim->beam_aperture) {
        jabs_message(MSG_ERROR, "Aperture could not be parsed.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc) {
        jabs_message(MSG_ERROR, "Unexpected extra arguments (%i), starting with %s.\n", argc, argv[0]);
        if(fit->sim->beam_aperture->type == APERTURE_NONE) {
            jabs_message(MSG_INFO, "Aperture type not defined. Allowed types:");
            for(const jibal_option *o = aperture_option; o->s; o++) {
                jabs_message(MSG_INFO, " %s", o->s);
            }
            jabs_message(MSG_INFO, "\n");
        } else {
            jabs_message(MSG_ERROR, "Aperture type was %s. Aperture removed.\n", aperture_name(fit->sim->beam_aperture));
        }
        aperture_free(fit->sim->beam_aperture);
        fit->sim->beam_aperture = NULL;
        return SCRIPT_COMMAND_FAILURE;
    }
    return argc_orig - argc;
}

script_command_status script_set_detector(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    size_t i_det = 0;
    const int argc_orig = argc;
    if(script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    s->i_det_active = i_det; /* This is not used beyond a single command currently, but this is possible in the future. */
    if(argc == 0) {
        return SCRIPT_COMMAND_NOT_FOUND;
    }
    return argc_orig - argc;
}

script_command_status script_set_detector_aperture(struct script_session *s, int argc, char *const *argv) {
    const int argc_orig = argc;
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Usage: set detector aperture <type> {width|height|diameter} ...\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    detector *det = sim_det(s->fit->sim, s->i_det_active);
    if(!det) {
        jabs_message(MSG_ERROR, "No detector(s)\n");
        return SCRIPT_COMMAND_SUCCESS;
    }
    if(detector_aperture_set_from_argv(s->jibal, det, &argc, &argv)) { /* TODO: handle parsing with subcommands */
        return SCRIPT_COMMAND_FAILURE;
    }
    return argc_orig - argc;
}

script_command_status script_set_detector_calibration(struct script_session *s, int argc, char *const *argv) {
    (void) s;
    (void) argc;
    (void) argv;
    const int argc_orig = argc;
    if(argc == 0) {
        return SCRIPT_COMMAND_NOT_FOUND;
    }
    const jibal_element *e = jibal_element_find(s->jibal->elements, argv[0]);
    if(e) {
        s->Z_active = e->Z;
        argc--;
    } else {
        s->Z_active = JIBAL_ANY_Z;
    }
    DEBUGMSG("Active Z is now %i", s->Z_active);
    return argc_orig - argc;
}

script_command_status script_set_detector_foil(struct script_session *s, int argc, char *const *argv) {
    const int argc_orig = argc;
    detector *det = sim_det(s->fit->sim, s->i_det_active);
    if(!det) {
        jabs_message(MSG_ERROR, "No detector(s)\n");
        return SCRIPT_COMMAND_SUCCESS;
    }
    if(detector_foil_set_from_argv(s->jibal, det, &argc, &argv)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    return argc_orig - argc;
}

script_command_status script_set_detector_calibration_poly(struct script_session *s, int argc, char *const *argv) {
    const int argc_orig = argc;
    detector *det = sim_det(s->fit->sim, s->i_det_active);
    if(!det) {
        jabs_message(MSG_ERROR, "No detector(s)\n");
        return SCRIPT_COMMAND_SUCCESS;
    }
    if(argc < 1) {
        jabs_message(MSG_ERROR,
                     "Usage: set detector calibration poly <n> <p_0> <p_1> ... <p_(n+1)>\nExample: set calibration poly 2 10keV 1.0keV 0.001keV\nThe example sets a second degree (quadratic 3 parameters) polynomial calibration.\nThe first parameter (p_0) is the constant term.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    size_t n = strtoull(argv[0], NULL, 10);
    if(n == 0) {
        jabs_message(MSG_ERROR, "Calibration should be a zero degree polynomial? Just no.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    argc--;
    argv++;
    if(argc < (int) (n + 1)) {
        jabs_message(MSG_ERROR, "Not enough parameters for a %zu degree polynomial. Expected %zu.\n", n, n + 1);
        return SCRIPT_COMMAND_FAILURE;
    }
    calibration *c = calibration_init_poly(n);
    calibration_copy_params(c, detector_get_calibration(det, s->Z_active)); /* This, de facto, only copies resolution from old calibration, since the rest are overwritten very soon. */
    for(int i = 0; i <= (int) n; i++) {
        if(jabs_unit_convert(s->jibal->units, JIBAL_UNIT_TYPE_ANY, argv[0], calibration_get_param_ref(c, i)) < 0) {
            calibration_free(c);
            return SCRIPT_COMMAND_FAILURE;
        }
        argc--;
        argv++;
    }
    calibration_params_poly *params = (calibration_params_poly *) c->params;
    for(int i = (int) n; i > 0; i--) { /* Reduce the degree of polynomial to match reality */
        if(calibration_get_param(c, i) == 0.0) {
            params->n--;
        } else {
            break;
        }
    }
    if(params->n == 0) {
        jabs_message(MSG_ERROR, "Calibration should be a zero degree polynomial? Just no.\n");
        calibration_free(c);
        return SCRIPT_COMMAND_FAILURE;
    }
    detector_set_calibration_Z(s->jibal->config, det, c, s->Z_active);
    return argc_orig - argc;
}

script_command_status script_set_detector_name(struct script_session *s, int argc, char * const *argv) {
    detector *det = sim_det(s->fit->sim, s->i_det_active);
    if(!det) {
        jabs_message(MSG_ERROR, "No detector(s)\n");
        return SCRIPT_COMMAND_SUCCESS;
    }
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Usage: set detector {<detector>} name <name>\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(detector_set_name(det, argv[0])) {
        jabs_message(MSG_ERROR, "Could not set name \"%s\". Please note some names are forbidden because they could be confused with detector numbers and commands.\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    return 1;
}

script_command_status script_set_sample(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 2) {
        jabs_message(MSG_ERROR, "Usage: set sample {nosimplify} <sample description>\nExample: set sample TiO2 1000tfu Si 10000tfu\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    const int argc_orig = argc;
    sample_model *sm_new = sample_model_from_argv(fit->jibal, &argc, &argv);
    if(!sm_new) {
        jabs_message(MSG_WARNING, "Setting sample fails due to parsing error.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    int argc_consumed = argc_orig - argc;
    if(fit_data_set_sample_model(fit, sm_new)) {
        jabs_message(MSG_WARNING, "Setting sample fails due to sample model sanity check.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(s->fit->sim->n_reactions > 0) {
        jabs_message(MSG_WARNING, "Reactions were reset automatically, since the sample was changed.\n");
    }
    return argc_consumed;
}

script_command_status script_set_stopping(struct script_session *s, int argc, char *const *argv) {
    const int argc_orig = argc;
    if(argc < 3) {
        jabs_message(MSG_ERROR, "Usage: set stopping <file> <element or Z1> <element or Z2> {<element or Z2>}\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(!s->jibal->gsto) {
        fprintf(stderr, "Unusual failure related to GSTO. Please report this.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    gsto_file_t *file = NULL;
    for(size_t i_file = 0; i_file < s->jibal->gsto->n_files; i_file++) {
        if(strcmp(s->jibal->gsto->files[i_file].name, argv[0]) == 0) {
            file = &(s->jibal->gsto->files[i_file]);
            break;
        }
    }
    if(!file) {
        jabs_message(MSG_ERROR, "Error: \"%s\" is not a file recognized by GSTO.\n", argv[0]);
        if(s->jibal->gsto->n_files) {
            jabs_message(MSG_INFO, "List of files:");
            for(size_t i_file = 0; i_file < s->jibal->gsto->n_files; i_file++) {
                jabs_message(MSG_INFO, " %s", s->jibal->gsto->files[i_file].name);
            }
            jabs_message(MSG_INFO, "\n");
        } else {
            jabs_message(MSG_INFO, "There are no GSTO files, they should be listed in this file: \"%s\"\n",
                         s->jibal->config->files_file);
        }
        return SCRIPT_COMMAND_FAILURE;
    }
    argc--;
    argv++;
    const jibal_element *first = jibal_element_find(s->jibal->elements, argv[0]);
    if(!first) {
        return SCRIPT_COMMAND_FAILURE;
    }
    argc--;
    argv++;
    while(argc) {
        const jibal_element *second = jibal_element_find(s->jibal->elements, argv[0]);
        if(!second)
            break;
        argc--;
        argv++;
        if(!jibal_gsto_assign(s->jibal->gsto, first->Z, second->Z, file)) {
            jabs_message(MSG_ERROR, "Could not assign stopping for Z1 = %i in Z2 = %i to file \"%s\".\n",
                         first->Z, second->Z, file->name);
            return SCRIPT_COMMAND_FAILURE;
        }
    }
    return argc_orig - argc;
}

script_command_status script_set_verbosity(struct script_session *s, int argc, char * const *argv) {
    (void) s;
    if(argc < 1) { /* argc_min set elsewhere, this is redundant, no error reporting */
        return EXIT_FAILURE;
    }
    size_t tmp = jabs_message_verbosity;
    if(jabs_str_to_size_t(argv[0], &tmp) < 0) {
        return EXIT_FAILURE;
    }
    if(tmp > MSG_ERROR) {
        jabs_message(MSG_ERROR, "The verbosity level %zu is too high, maximum is %zu.\n", tmp, MSG_ERROR);
        return EXIT_FAILURE;
    }
    jabs_message_verbosity = tmp;
    jabs_message(MSG_INFO, "Verbosity level set to \"%s\".\n", jabs_message_level_str(jabs_message_verbosity));
    return 1;
}

script_command_status script_set_reaction_yield(struct script_session *s, int argc, char * const *argv) {
    (void) s;
    if(argc < 1) { /* argc_min set elsewhere, this is redundant, no error reporting */
        return EXIT_FAILURE;
    }
    size_t tmp = jabs_message_verbosity;
    if(jabs_str_to_size_t(argv[0], &tmp) < 0) {
        return EXIT_FAILURE;
    }
    if(tmp > MSG_ERROR) {
        jabs_message(MSG_ERROR, "The verbosity level %zu is too high, maximum is %zu.\n", tmp, MSG_ERROR);
        return EXIT_FAILURE;
    }
    jabs_message_verbosity = tmp;
    jabs_message(MSG_INFO, "Verbosity level set to \"%s\".\n", jabs_message_level_str(jabs_message_verbosity));
    return 1;
}

script_command_status script_test_reference(struct script_session *s, int argc, char *const *argv) {
    const struct fit_data *fit = s->fit;
    const int argc_orig = argc;
    size_t i_det = 0;
    if(!fit->ref) {
        jabs_message(MSG_ERROR, "No reference spectrum loaded. Use 'load reference' to do so.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc < 2) {
        jabs_message(MSG_ERROR, "Usage: test reference {detector number} <range> <tolerance>\nTests if simulated spectrum is within tolerance to a reference spectrum.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    const jabs_histogram *sim = fit_data_histo_sum(fit, i_det);
    if(!sim) {
        jabs_message(MSG_ERROR, "No simulation.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    roi r = {.i_det = i_det};
    if(fit_set_roi_from_string(&r, argv[0])) {
        jabs_message(MSG_ERROR, "Could not parse range.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    double tolerance = strtod(argv[1], NULL);
    double error;
    argc -= 2;
    argv += 2;
    int return_value = argc_orig - argc;
    if(jabs_histogram_compare(sim, fit->ref, r.low, r.high, &error)) { /* Failed, test fails */
        jabs_message(MSG_ERROR, "Test of detector %zu from %zu to %zu failed. Is range valid?\n", r.i_det + 1, r.low, r.high);
        return_value = SCRIPT_COMMAND_FAILURE;
    } else {
        jabs_message(MSG_INFO, "Test of simulated spectrum to reference from %zu to %zu. Error %e.\n", r.low, r.high, error);
        if(error > tolerance) {
            jabs_message(MSG_ERROR, "Test failed.\n");
            return_value = SCRIPT_COMMAND_FAILURE;
        } else {
            jabs_message(MSG_INFO, "Test passed.\n");
        }
    }
    return return_value;
}

script_command_status script_test_roi(struct script_session *s, int argc, char *const *argv) {
    const struct fit_data *fit = s->fit;
    const int argc_orig = argc;
    size_t i_det = 0;
    int exp = FALSE;

    if(script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc && strncmp(argv[0], "exp", 3) == 0) { /* "exp" argument given, compare against experimental spectrum */
        exp = TRUE;
        argc--;
        argv++;
    }

    if(argc < 3) {
        jabs_message(MSG_ERROR, "Usage: test roi {<detector>} {exp} <range> <sum> <tolerance>\nTests if simulated (or experimental) ROI sum is within relative tolerance.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    roi r = {.i_det = i_det};
    if(fit_set_roi_from_string(&r, argv[0])) {
        jabs_message(MSG_ERROR, "Could not parse range.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    double sum_ref = strtod(argv[1], NULL);
    double tolerance = strtod(argv[2], NULL);
    jabs_histogram *h = exp ? fit_data_exp(fit, i_det) : fit_data_histo_sum(fit, i_det);
    if(!h) {
        jabs_message(MSG_ERROR, "No histogram.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    double sum = jabs_histogram_roi(h, r.low, r.high);
    argc -= 3;
    argv += 3;
    double rel_err = fabs(1.0 - sum / sum_ref);
    jabs_message(MSG_INFO, "Test of ROI from %zu to %zu. Sum %.10g. Relative error %e.\n", r.low, r.high, sum, rel_err);
    if(rel_err < tolerance) {
        jabs_message(MSG_INFO, "Test passed.\n");
    } else {
        jabs_message(MSG_ERROR, "Test failed.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    return argc_orig - argc;
}

script_command_status script_split_sample_elements(struct script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    struct fit_data *fit = s->fit;
    sample_model *sm = fit->sm;
    fit->sm = sample_model_split_elements(sm);
    sample_model_free(sm);
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_add_reaction(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 2) {
        jabs_message(MSG_ERROR, "Usage: add reaction <type> <isotope> {cs <cross section model>} {min <min energy>} {max <max energy}.\n Type should be one of the following: RBS RBS- ERD\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    const int argc_orig = argc;
    reaction *r = sim_reaction_make_from_argv(fit->jibal, fit->sim, &argc, &argv);
    int argc_consumed = argc_orig - argc;
    if(!r) {
        jabs_message(MSG_ERROR, "Could not make a reaction based on given description.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(sim_reactions_add_reaction(fit->sim, r, FALSE)) {
        return SCRIPT_COMMAND_FAILURE;
    } else {
        return argc_consumed;
    }
}

script_command_status script_add_reactions(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    if(!fit->sm) {
        jabs_message(MSG_ERROR,
                     "Cannot add reactions before sample has been set (I need to know which reactions to add!).\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc > 0) {
        reaction_type type = reaction_type_from_string(argv[0]);
        if(type == REACTION_NONE) {
            jabs_message(MSG_ERROR, "Could not determine what reaction type %s is? Did you mean \"RBS\"?\n", argv[0]);
            return SCRIPT_COMMAND_FAILURE;
        }
        if(type == REACTION_RBS || type == REACTION_RBS_ALT || type == REACTION_ERD) {
            sim_reactions_add_auto(fit->sim, fit->sm, type, sim_cs(fit->sim, type), FALSE);
        } else {
            jabs_message(MSG_ERROR, "Adding reactions of type %s is either not implemented or it is done by another command.\n", reaction_type_to_string(type));
            return SCRIPT_COMMAND_FAILURE;
        }
        return 1;
    }
    if(fit->sim->rbs) {
        sim_reactions_add_auto(fit->sim, fit->sm, REACTION_RBS,
                               sim_cs(fit->sim, REACTION_RBS), FALSE); /* TODO: loop over all detectors and add reactions that are possible (one reaction for all detectors) */
        sim_reactions_add_auto(fit->sim, fit->sm, REACTION_RBS_ALT,
                               sim_cs(fit->sim, REACTION_RBS_ALT), FALSE);
    }

    if(sim_do_we_need_erd(fit->sim)) {
        sim_reactions_add_auto(fit->sim, fit->sm, REACTION_ERD, sim_cs(fit->sim, REACTION_ERD), FALSE);
    }
    return 0;
}

script_command_status script_add_detector_default(script_session *s, int argc, char *const *argv) {
    (void) argc; /* Doesn't consume arguments */
    (void) argv;
    struct fit_data *fit = s->fit;
    if(fit_data_add_det(fit, detector_default(NULL))) {
        return SCRIPT_COMMAND_FAILURE;
    } else {
        return 0;
    }
}

script_command_status script_add_fit_range(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    size_t i_det = 0;
    const int argc_orig = argc;
    script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det);
    int n_ranges = 0;
    if(argc < 1) {
        jabs_message(MSG_ERROR,
                     "Usage: add fit_range {detector} <range> {<range> <range> ...}\nExample: add fit_range 1 [400:900] [980:1200]\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    while(argc > 0) {
        roi r = {.i_det = i_det};
        if(fit_set_roi_from_string(&r, argv[0])) {
            if(n_ranges == 0) {
                jabs_message(MSG_ERROR, "No ranges added! Failed on parsing \"%s\".\n", argv[0]);
                return SCRIPT_COMMAND_FAILURE;
            } else {
                break;
            }
        }
        if(fit_data_fit_range_add(fit, &r)) {
            jabs_message(MSG_ERROR, "Could not add range \"%s\" for some reason.\n", argv[0]);
            return SCRIPT_COMMAND_FAILURE;
        }
        n_ranges++;
        argc--;
        argv++;
    }
    return argc_orig - argc;
}

script_command_status script_help(script_session *s, int argc, char *const *argv) {
    (void) s;
    static const struct help_topic topics[] = {
            {"help", "This is help on help.\nHelp is available on following topics:\n"},
            {NULL, NULL}
    };
    if(argc == 0) {
        jabs_message(MSG_INFO, "Type help [topic] for information on a particular topic or \"help help\" for help on help.\n\n");

        return SCRIPT_COMMAND_NOT_FOUND;
    }

    int found = 0;
    for(const struct help_topic *t = topics; t->name != NULL; t++) {
        if(strcmp(t->name, argv[0]) == 0) {
            found++;
            jabs_message(MSG_INFO, "%s", t->help_text);
            if(strcmp(t->name, "help") == 0) {
                size_t i = 0;
                for(const struct help_topic *t2 = topics; t2->name != NULL; t2++) {
                    i++;
                    jabs_message(MSG_INFO, "%18s", t2->name);
                    if(i % 4 == 0) {
                        jabs_message(MSG_INFO, "\n");
                    }
                }
                jabs_message(MSG_INFO, "\nTry also \"help commands\" for a list of possible commands. Try running a command without arguments to get a brief help on usage.");
            }
            jabs_message(MSG_INFO, "\n");
            break;
        }
    }

    if(!found) {
        return SCRIPT_COMMAND_NOT_FOUND;
    }
    return 1; /* TODO: this stuff only works with one argument, right? */
}

script_command_status script_help_version(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    (void) s;
    if(git_populated()) {
        jabs_message(MSG_INFO, "JaBS version %s (git: %s).\n", jabs_version_simple(), jabs_version());
        jabs_message(MSG_INFO, "This version of JaBS is compiled from a git repository (branch %s%s).\n", git_branch(), git_dirty() ? ", dirty" : "");
        jabs_message(MSG_INFO, "Git commit %s dated %s.\n", git_commit_sha1(), git_commit_date());
    } else {
        jabs_message(MSG_INFO, "JaBS version %s (not from git repository).\n", jabs_version());
    }
#ifdef __DATE__
    jabs_message(MSG_INFO,  "Compiled on %s\n", __DATE__);
#endif
#ifdef __VERSION__
    jabs_message(MSG_INFO,  "Compiled with compiler version: %s\n", __VERSION__);
#endif
#ifdef _MSC_VER
    jabs_message(MSG_INFO, "Compiled with MSVC version %i\n", _MSC_VER);
#endif
#ifdef JABS_PLUGINS
    jabs_message(MSG_INFO,  "Plugin support enabled.\n");
#endif
#ifdef _OPENMP
    jabs_message(MSG_INFO,  "OpenMP: %d\n", _OPENMP);
#endif
    jabs_message(MSG_INFO, "GSL: %s (compile time: %s)\n", gsl_version, GSL_VERSION);
    jabs_message(MSG_INFO, "JIBAL: %s (compile time: %s)\n", jibal_version(), JIBAL_VERSION);
    jabs_message(MSG_INFO, "libxml2 compile time: %s\n", LIBXML_DOTTED_VERSION);
    LIBXML_TEST_VERSION
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_help_commands(script_session *s, int argc, char *const *argv) {
    (void) argv;
    if(argc == 0) {
        jabs_message(MSG_INFO, "The following commands are available:\n");
        script_print_command_tree(s->commands);
        return SCRIPT_COMMAND_SUCCESS;
    }
    return SCRIPT_COMMAND_NOT_FOUND;
}

#ifdef JABS_PLUGINS
script_command_status script_identify_plugin(struct script_session *s, int argc, char * const *argv) {
    (void) s;
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Usage: identify plugin <path>\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    jabs_plugin *plugin = jabs_plugin_open(argv[0]);
    if(!plugin) {
        jabs_message(MSG_ERROR, "Plugin %s not found or could not be opened.\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    jabs_message(MSG_INFO, "Plugin identifies as \"%s\" version \"%s\", type of plugin is %s.\n", plugin->name, plugin->version, jabs_plugin_type_string(plugin->type));
    jabs_plugin_close(plugin);
    return 1;
}

script_command_status script_load_reaction_plugin(script_session *s, int argc, char *const *argv) {
    const int argc_orig = argc;
    struct fit_data *fit = s->fit;
    if(argc < 2) {
        jabs_message(MSG_ERROR, "Usage: load reaction plugin <file> <target> ...\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    const jibal_isotope *target = jibal_isotope_find(s->jibal->isotopes, argv[1], 0, 0);
    if(!target) {
        jabs_message(MSG_ERROR, "Could not find isotope matching \"%s\".\n", argv[1]);
        return SCRIPT_COMMAND_FAILURE;
    }
    const char *filename = argv[0];
    jabs_plugin *plugin = jabs_plugin_open(filename);
    if(!plugin) {
        jabs_message(MSG_ERROR, "Could not load plugin from file \"%s\".\n", filename);
        return EXIT_FAILURE;
    }
    reaction *r = reaction_make(fit->sim->beam_isotope, target, REACTION_PLUGIN, JABS_CS_NONE);
    if(!r) {
        jabs_message(MSG_ERROR, "Could not make a new reaction.\n");
        jabs_plugin_close(plugin);
        return SCRIPT_COMMAND_FAILURE;
    }
    r->plugin = plugin;
    argc -= 2;
    argv += 2;
    jabs_plugin_reaction *pr = jabs_plugin_reaction_init(plugin, s->jibal->isotopes, fit->sim->beam_isotope, target, &argc, &argv);
    if(!pr) {
        jabs_message(MSG_ERROR, "Plugin failed to initialize.\n");
        reaction_free(r);
        return SCRIPT_COMMAND_FAILURE;

    }
    r->plugin_r = pr;
    r->product = pr->product;
    r->residual = pr->product_heavy;
    r->E_min = pr->E_min;
    r->E_max = pr->E_max;
    r->filename = strdup(plugin->filename);
    r->Q = pr->Q;
    r->yield = pr->yield;
    reaction_generate_name(r);
    sim_reactions_add_reaction(fit->sim, r, FALSE);
    return argc_orig; /* TODO: always consumes all arguments */
}
#endif

script_command_status script_cwd(struct script_session *s, int argc, char *const *argv) {
    (void) s;
    (void) argv;
    if(argc != 0) {
        return SCRIPT_COMMAND_NOT_FOUND;
    }
    char *cwd = getcwd(NULL, 0);
    jabs_message(MSG_INFO, "%s\n", cwd);
    free(cwd);
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_cd(struct script_session *s, int argc, char *const *argv) {
    (void) s;
    const int argc_orig = argc;
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Usage: cd <path>\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(chdir(argv[0])) {
        jabs_message(MSG_ERROR, "Could not change directory to \"%s\"\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    argc--;
    argv++;
    return argc_orig - argc;
}

script_command_status script_idf2jbs(struct script_session *s, int argc, char * const *argv) {
    const int argc_orig = argc;
    (void) s;
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Usage: idf2jbs <idf file>\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    char *filename_out = NULL;
    idf_error idferr = idf_parse(argv[0], &filename_out);
    if(idferr == IDF2JBS_SUCCESS) {
        jabs_message(MSG_INFO, "Success. Wrote script to file \"%s\"\n", filename_out);
    } else {
        jabs_message(MSG_ERROR, "IDF2JBS failed with error code %i (%s).\n", idferr, idf_error_code_to_str(idferr));
    }
    free(filename_out);
    argc--;
    argv++;
    return argc_orig - argc;
}

script_command_status script_kinematics(struct script_session *s, int argc, char * const *argv) {
    const int argc_orig = argc;
    struct fit_data *fit = s->fit;
    const simulation *sim = s->fit->sim;
    if(argc < 1) {
        jabs_message(MSG_ERROR, "Usage: kinematics <type> <target atom>\nSyntax of \"add reaction\" is used.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(sim->n_det == 0) {
        jabs_message(MSG_ERROR, "No detectors defined.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    reaction *r = sim_reaction_make_from_argv(fit->jibal, fit->sim, &argc, &argv);
    if(!r) {
        jabs_message(MSG_ERROR, "Reaction could not be initialized.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    jabs_message(MSG_INFO, "Detector name | theta |  product |   energy | channel | cross section\n");
    jabs_message(MSG_INFO, "              |   deg |          |      keV |         |         mb/sr\n");
    for(size_t i_det = 0; i_det < sim->n_det; i_det++) {

        const detector *det = sim->det[i_det];
        if(reaction_is_possible(r, sim->params, det->theta)) {
            sim_reaction *sim_r = calloc(1, sizeof(sim_reaction)); /* We skip initializing the complete reaction, just use parts of it */
            sim_r->r = r;
            sim_reaction_set_cross_section_by_type(sim_r);
            sim_r->emax_incident = sim->beam_E;
            sim_reaction_recalculate_internal_variables(sim_r, sim->params, det->theta, sim->emin, sim->beam_E, sim->beam_E);
            sim_reaction_recalculate_screening_table(sim_r); /* TODO: only when necessary (Universal) */
            double E = reaction_product_energy(r, det->theta, sim->beam_E);
            const jibal_isotope *product = r->product;
            if(!product) {
                continue;
            }
            size_t ch = calibration_inverse(detector_get_calibration(det, product->Z), E, CHANNELS_ABSOLUTE_MAX);
            double sigma = sim_r->cross_section ? sim_r->cross_section(sim_r, sim->beam_E) : 0.0;
            jabs_message(MSG_INFO, "%13s | %5g | %8s | %8g | %7zu | %13g\n",
                         detector_name(det),
                         det->theta / C_DEG,
                         product->name,
                         E / C_KEV,
                         ch,
                         sigma / C_MB_SR);
            free(sim_r);
        } else {

        }
    }
    reaction_free(r);
    return argc_orig - argc;
}
