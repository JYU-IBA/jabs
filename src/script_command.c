/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
#include <direct.h> // getcwd
#else
#include <unistd.h>
#endif
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


int script_prepare_sim_or_fit(script_session *s) {
    fit_data *fit = s->fit;
    if(!fit->sm) {
        jabs_message(MSG_ERROR, stderr, "No sample has been defined!\n");
        return -1;
    }
    if(!fit->sim->beam_isotope) {
        jabs_message(MSG_ERROR, stderr, "No ion has been defined!\n");
        return -1;
    }
    if(!fit->sim->det || fit->sim->n_det == 0) {
        jabs_message(MSG_ERROR, stderr, "No detector has been defined!\n");
        return -1;
    }
    if(sim_sanity_check(fit->sim)) {
        jabs_message(MSG_ERROR, stderr, "Simulation failed sanity check.\n");
        return -1;
    }
    fit_data_workspaces_free(s->fit);
    sample_free(fit->sim->sample);
#ifdef DEBUG
    fprintf(stderr, "Original sample model:\n");
    sample_model_print(NULL, fit->sm);
#endif
    fit->sim->sample = sample_from_sample_model(fit->sm);
    if(!fit->sim->sample) {
        jabs_message(MSG_ERROR, stderr,
                     "Could not make a sample based on model description. This should never happen.\n");
        return -1;
    }
    if(fit->sim->n_reactions == 0) {
        jabs_message(MSG_WARNING, stderr,
                     "No reactions, adding some automatically. Please be aware there are commands called \"reset reactions\" and \"add reactions\".\n");
        if(fit->sim->rbs) {
            sim_reactions_add_auto(fit->sim, fit->sm, REACTION_RBS, sim_cs(fit->sim, REACTION_RBS));
            sim_reactions_add_auto(fit->sim, fit->sm, REACTION_RBS_ALT, sim_cs(fit->sim, REACTION_RBS_ALT));
            /* TODO: loop over all detectors and add reactions that are possible (one reaction for all detectors) */
        }
        if(sim_do_we_need_erd(fit->sim)) {
            sim_reactions_add_auto(fit->sim, fit->sm, REACTION_ERD, sim_cs(fit->sim, REACTION_ERD));
        }
    }
    if(fit->sim->n_reactions == 0) {
        jabs_message(MSG_ERROR, stderr, "No reactions. Nothing to do.\n");
        return EXIT_FAILURE;
    }
    sim_sort_reactions(fit->sim);
    jabs_message(MSG_INFO, stderr, "Simplified sample model for simulation:\n");
    sample_print(NULL, fit->sim->sample, TRUE);

    reactions_print(stderr, fit->sim->reactions, fit->sim->n_reactions);

    if(assign_stopping(fit->jibal->gsto, fit->sim)) {
        jabs_message(MSG_ERROR, stderr,
                     "Could not assign stopping or straggling. Failure. Provide more data, check that JIBAL Z2_max is sufficiently large (currently %i) or disable unwanted reactions (e.g. ERD).\n",
                     s->jibal->config->Z_max);
        return -1;
    }
    script_show_stopping(s, 0, NULL);
    jibal_gsto_print_files(fit->jibal->gsto, TRUE); /* TODO: this don't use jabs_message() */
    jabs_message(MSG_VERBOSE, stderr, "Loading stopping data.\n");
    jibal_gsto_load_all(fit->jibal->gsto);
#ifdef DEBUG
    fprintf(stderr, "Updating calculation params before sim/fit\n");
#endif
    sim_calc_params_update(fit->sim->params);
    sim_print(fit->sim);
    s->start = clock();
    return 0;
}

int script_finish_sim_or_fit(script_session *s) {
    s->end = clock();
    double cputime_total = (((double) (s->end - s->start)) / CLOCKS_PER_SEC);
    jabs_message(MSG_INFO, stderr, "\n...finished! Total CPU time: %.3lf s.\n", cputime_total);
#ifdef CLEAR_GSTO_ASSIGNMENTS_WHEN_FINISHED
    jibal_gsto_assign_clear_all(s->fit->jibal->gsto); /* Is it necessary? No. Here? No. Does it clear old stuff? Yes. */
#endif

    struct fit_data *fit = s->fit;

    if(fit->sim->n_det == 1) { /* TODO: This is primarily used for command line mode, but multidetector mode could be supported. */
        size_t i_det = 0;
        sim_workspace *ws = fit_data_ws(fit, i_det);
        if(ws) {
            if(s->output_filename) {
                if(print_spectra(s->output_filename, ws, fit_data_exp(fit, i_det))) {
                    jabs_message(MSG_ERROR, stderr, "Could not save spectra of detector %zu to file \"%s\"\n", i_det,
                                 s->output_filename);
                    return EXIT_FAILURE;
                }
            }
        }
    }
    return 0;
}


void script_command_not_found(const char *cmd, const script_command *c_parent) {
    if(c_parent) {
        if(cmd) {
            jabs_message(MSG_ERROR, stderr, "Sub-command \"%s\" is invalid!\n\n", cmd);
        } else {
            jabs_message(MSG_ERROR, stderr, "Not enough arguments!\n\n");
        }
        if(c_parent->subcommands) {
            size_t matches = script_command_print_possible_matches_if_ambiguous(c_parent->subcommands, cmd);
            if(matches == 0) {
                jabs_message(MSG_ERROR, stderr, "Following subcommands of \"%s\" are recognized:\n", c_parent->name);
                script_commands_print(stderr, c_parent->subcommands);
            }
        }
    } else {
        if(cmd) {
            jabs_message(MSG_ERROR, stderr, "Invalid command or argument: \"%s\".\n", cmd);
        } else {
            jabs_message(MSG_ERROR, stderr, "What?\n", cmd);
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
    if(fit_data_workspaces_init(fit)) {
        jabs_message(MSG_ERROR, stderr, "Could not initialize simulation workspace(s).\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    sim_calc_params_print(fit->sim->params);
    for(size_t i_det = 0; i_det < fit->sim->n_det; i_det++) {
        if(simulate_with_ds(fit->ws[i_det])) {
            jabs_message(MSG_ERROR, stderr, "Simulation failed.\n");
            return SCRIPT_COMMAND_FAILURE;
        }
    }
    script_finish_sim_or_fit(s);
    return argc_orig - argc;
}

int foo(fit_params *params, const char *fit_vars) {
#ifdef DEBUG
    fprintf(stderr, "fitvars = %s\n", fit_vars);
#endif
    int status = EXIT_SUCCESS;
    if(!fit_vars)
        return EXIT_FAILURE;
    char *token, *s, *s_orig;
    s_orig = s = strdup(fit_vars);
    while((token = strsep_with_quotes(&s, ",")) != NULL) { /* parse comma separated list of parameters to fit */
        if(fit_params_enable(params, token, TRUE) == 0) {
            jabs_message(MSG_ERROR, stderr, "No matches for %s. See 'show fit variables' for a list of possible fit variables.\n", token);
            status = EXIT_FAILURE;
        }
        if(status == EXIT_FAILURE)
            break;
    }
    free(s_orig);
    return status;
}

script_command_status script_fit(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit_data = s->fit;
    if(argc != 1) {
        fprintf(stderr, "Usage: fit [fitvar1,fitvar2,...]\nSee 'show fit variables' for a list of possible fit variables.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    fit_params_free(fit_data->fit_params);

    fit_params *p_all = fit_params_all(fit_data);
    if(foo(p_all, argv[0])) {
        jabs_message(MSG_ERROR, stderr, "Error in adding fit parameters.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    fit_params_update(p_all);
    fit_params_print(p_all, FALSE, NULL);
    fit_data->fit_params = p_all;

    if(fit_data->fit_params->n_active == 0) {
        jabs_message(MSG_ERROR, stderr, "No parameters for fit.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    jabs_message(MSG_INFO, stderr, "%zu fit parameters, %zu active.\n", fit_data->fit_params->n, fit_data->fit_params->n_active);
    if(!fit_data->exp) { /* TODO: not enough to check this */
        jabs_message(MSG_ERROR, stderr, "No experimental spectrum set.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(fit_data->n_fit_ranges == 0) {
        jabs_message(MSG_ERROR, stderr, "No fit range(s) given.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(script_prepare_sim_or_fit(s)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    fit_data->fit_iter_callback = s->fit_iter_callback;
    if(fit(fit_data) < 0) {
        return SCRIPT_COMMAND_FAILURE;
    }
    script_finish_sim_or_fit(s);
#ifdef PRINT_SIM_AFTER_FIT
    jabs_message(MSG_INFO, stderr, "\nFinal parameters:\n");
    simulation_print(stderr, fit_data->sim);
#endif
    jabs_message(MSG_INFO, stderr, "\nFinal profile:\n");
    sample_print(NULL, fit_data->sim->sample, FALSE);
    sample_areal_densities_print(stderr, fit_data->sim->sample, FALSE);
    jabs_message(MSG_INFO, stderr, "\nFinal layer thicknesses:\n");
    sample_print_thicknesses(NULL, fit_data->sim->sample);
    jabs_message(MSG_INFO, stderr, "\nFinal sample model:\n");
    sample_model_print(NULL, fit_data->sm);
    jabs_message(MSG_INFO, stderr, "\n");
    fit_stats_print(stderr, &fit_data->stats);
    return 1;
}

script_command_status script_save_bricks(script_session *s, int argc, char *const *argv) {
    size_t i_det = 0;
    const int argc_orig = argc;
    struct fit_data *fit = s->fit;
    if(script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det) || argc != 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: save bricks [detector] file\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(print_bricks(argv[0], fit_data_ws(fit, i_det))) {
        jabs_message(MSG_ERROR, stderr,
                     "Could not save bricks of detector %zu to file \"%s\"! There should be %zu detector(s).\n",
                     i_det + 1, argv[0], fit->sim->n_det);
        return SCRIPT_COMMAND_FAILURE;
    }
    argc--;
    argv++;
    return argc_orig - argc;
}

script_command_status script_save_spectra(script_session *s, int argc, char *const *argv) {
    size_t i_det = 0;
    struct fit_data *fit = s->fit;
    const int argc_orig = argc;
    if(script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det) || argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: save spectra {<detector>} file\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Not enough arguments for save spectra.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(print_spectra(argv[0], fit_data_ws(fit, i_det), fit_data_exp(fit, i_det))) {
        jabs_message(MSG_ERROR, stderr,
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
        jabs_message(MSG_ERROR, stderr, "Usage: save sample <file>\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(!fit_data->sm) {
        jabs_message(MSG_ERROR, stderr, "No sample set.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(sample_model_print(argv[0], fit_data->sm)) {
        jabs_message(MSG_ERROR, stderr, "Could not write sample to file \"%s\".\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    return 1;
}

script_command_status script_save_calibrations(script_session *s, int argc, char *const *argv) {
    const int argc_orig = argc;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: save calibrations <file>\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    FILE *f = fopen_file_or_stream(argv[0], "w");
    if(!f) {
        jabs_message(MSG_ERROR, stderr, "Can not open file \"%s\" for writing.\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    argc--;
    for(size_t i = 0; i < s->fit->sim->n_det; i++) {
        const detector *det = sim_det(s->fit->sim, i);
        char *calib_str = calibration_to_string(det->calibration);
        char *reso_str = detector_resolution_to_string(det, JIBAL_ANY_Z);
        jabs_message(MSG_INFO, f, "set detector %zu calibration %s resolution %s\n", i + 1, calib_str, reso_str);
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
        jabs_message(MSG_ERROR, stderr, remove_reaction_usage);
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
        jabs_message(MSG_ERROR, stderr, remove_reaction_usage);
        return SCRIPT_COMMAND_FAILURE;
    }
    reaction_type type = reaction_type_from_string(argv[0]);
    const jibal_isotope *target = jibal_isotope_find(fit->jibal->isotopes, argv[1], 0, 0);
    if(type == REACTION_NONE) {
        jabs_message(MSG_ERROR, stderr, "This is not a valid reaction type: \"%s\".\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    if(!target) {
        jabs_message(MSG_ERROR, stderr, "This is not a valid isotope: \"%s\".\n", argv[1]);
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
    jabs_message(MSG_ERROR, stderr, "No matching reaction found!\n");
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
        jabs_message(MSG_ERROR, stderr, "Usage: roi <range> {<range> <range> ...}\nExample: roi [400:900] [980:1200]\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    while(argc > 0) {
        roi r = {.i_det = i_det};
        if(fit_set_roi_from_string(&r, argv[0])) {
            break;
        }
        fit_data_roi_print(stderr, s->fit, &r);
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
    int found = 0;
    const script_command *c_found = NULL;
    for(const script_command *c = commands; c; c = c->next) {
        if(strncmp(c->name, cmd_string, strlen(cmd_string)) == 0) {
            found++;
            c_found = c;
#ifdef DEBUG
            fprintf(stderr, "Candidate for \"%s\": \"%s\".\n", cmd_string, c->name);
#endif
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
        jabs_message(MSG_ERROR, stderr, "\"%s\" is ambiguous (%i matches):", cmd_string, found);
        for(const script_command *c = commands; c; c = c->next) {
            if(strncmp(c->name, cmd_string, strlen(cmd_string)) == 0) {
                jabs_message(MSG_ERROR, stderr, " %s", c->name);
            }
        }
        jabs_message(MSG_ERROR, stderr, "\n");
    }
    return found;
}

script_command_status script_set_var(struct script_session *s, jibal_config_var *var, int argc, char *const *argv) {
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Not enough arguments to set variable.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    jibal_config_var_set(s->jibal->units, var, argv[0], NULL);
    return 1; /* Number of arguments */
}

script_command_status script_set_detector_val(struct script_session *s, int val, int argc, char *const *argv) {
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Not enough arguments to set detector variable.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
#ifdef DEBUG
    fprintf(stderr, "Active detector is %zu.\n", s->i_det_active);
#endif
    detector *det = sim_det(s->fit->sim, s->i_det_active);
    if(!det) {
        jabs_message(MSG_ERROR, stderr, "Detector %zu does not exist.\n", s->i_det_active);
        return SCRIPT_COMMAND_FAILURE;
    }
    double *value_double = NULL;
    size_t *value_size = NULL;
    switch(val) {
        case 'b': /* beta */
            value_double = &(det->beta);
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
            break;
        case 't': /* type */
            det->type = jibal_option_get_value(detector_option, argv[0]);
            if(det->type == 0) {
                jabs_message(MSG_ERROR, stderr, "Detector type \"%s\" is none or unknown.\n", argv[0]);
                return SCRIPT_COMMAND_FAILURE;
            }
            break;
        case 'S': /* slope, this is for backwards compatibility (and ease of use with linear calibration) */
            value_double = calibration_get_param_ref(det->calibration, CALIBRATION_PARAM_SLOPE);
            break;
        case 'O': /* offset, this is for backwards compatibility */
            value_double = calibration_get_param_ref(det->calibration, CALIBRATION_PARAM_OFFSET);
            break;
        case 'r': /* resolution */
            value_double = calibration_get_param_ref(det->calibration, CALIBRATION_PARAM_RESOLUTION);
            break;
        case 's': /* solid */
            value_double = &(det->solid);
            break;
        case 'T': /* theta */
            value_double = &(det->theta);
            break;
        case 'l': /* length */
            value_double = &(det->length);
            break;
        case 'p': /* phi */
            value_double = &(det->phi);
            break;
        default:
            jabs_message(MSG_ERROR, stderr, "Unhandled value %i in script_set_detector_val. Report to developer.\n", val);
            return SCRIPT_COMMAND_FAILURE;
            break;
    }
    if(value_double) {
        *value_double = jibal_get_val(s->jibal->units, 0, argv[0]);
    } else if(value_size) {
        *value_size = strtoull(argv[0], NULL, 10);
    }
    return 1; /* Number of arguments */
}

script_command_status script_set_detector_calibration_val(struct script_session *s, int val, int argc, char *const *argv) {
    (void) argv;
#ifdef DEBUG
    fprintf(stderr, "Active detector is %zu.\n", s->i_det_active);
#endif
    detector *det = sim_det(s->fit->sim, s->i_det_active);
    if(!det) {
        jabs_message(MSG_ERROR, stderr, "Detector %zu does not exist.\n", s->i_det_active);
        return SCRIPT_COMMAND_FAILURE;
    }
    calibration *c;
    switch(val) { /* Handle cases where we don't expect (consume) arguments */
        case 'L': /* linear */
            c = calibration_init_linear();
            calibration_copy_params(c, detector_get_calibration(det, s->Z_active)); /* Copy parameters from old calibration, as much as possible */
            if(detector_set_calibration_Z(s->jibal->config, det, c, s->Z_active)) {
                jabs_message(MSG_ERROR, stderr, "Could not set linear calibration (element = %s).\n",
                             jibal_element_name(s->jibal->elements, s->Z_active));
                return SCRIPT_COMMAND_FAILURE;
            }
            return 0;
        default:
            break;
    }
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Not enough parameters to set detector calibration values.\n", val);
        return SCRIPT_COMMAND_FAILURE;
    }
    c = detector_get_calibration(det, s->Z_active);
    if(!c) {
        jabs_message(MSG_ERROR, stderr, "No calibration set for element = %s\n", jibal_element_name(s->jibal->elements, s->Z_active));
        return SCRIPT_COMMAND_FAILURE;
    }
    double value_double = jibal_get_val(s->jibal->units, 0, argv[0]);
    switch(val) {
        case 's': /* slope */
            if(calibration_set_param(c, CALIBRATION_PARAM_SLOPE, value_double)) {
                jabs_message(MSG_ERROR, stderr, "Can not set calibration slope.\n");
                return SCRIPT_COMMAND_FAILURE;
            }
            return 1;
        case 'o': /* offset */
            if(calibration_set_param(c, CALIBRATION_PARAM_OFFSET, value_double)) {
                jabs_message(MSG_ERROR, stderr, "Can not set calibration offset.\n");
                return SCRIPT_COMMAND_FAILURE;
            }
            return 1;
        case 'r': /* resolution */
            if(calibration_set_param(c, CALIBRATION_PARAM_RESOLUTION, value_double)) {
                jabs_message(MSG_ERROR, stderr, "Can not set calibration resolution.\n");
                return SCRIPT_COMMAND_FAILURE;
            }
            return 1;
        default:
            break;
    }
    jabs_message(MSG_ERROR, stderr, "Unhandled value %i in script_set_detector_calibration_val. Report to developer.\n", val);
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
        jabs_message(MSG_ERROR, stderr, "Not enough parameters to set fit values.\n");
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
        jabs_message(MSG_ERROR, stderr, "Not enough parameters to set simulation values.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    return SCRIPT_COMMAND_SUCCESS;
}


script_command_status script_show_var(struct script_session *s, jibal_config_var *var, int argc, char *const *argv) {
    (void) argv;
    (void) argc;
    (void) s;
    if(var->variable == NULL)
        return SCRIPT_COMMAND_FAILURE;
    switch(var->type) {
        case JIBAL_CONFIG_VAR_NONE:
            break;
        case JIBAL_CONFIG_VAR_PATH:
        case JIBAL_CONFIG_VAR_STRING:
            if(*((void **) var->variable) == NULL)
                break;
            jabs_message(MSG_INFO, stderr, "%s = %s\n", var->name, *((char **) var->variable));
            break;
        case JIBAL_CONFIG_VAR_BOOL:
            jabs_message(MSG_INFO, stderr, "%s = %s\n", var->name, *((int *) var->variable) ? "true" : "false");
            break;
        case JIBAL_CONFIG_VAR_INT:
            jabs_message(MSG_INFO, stderr, "%s = %i\n", var->name, *((int *) var->variable));
            break;
        case JIBAL_CONFIG_VAR_DOUBLE:
            jabs_message(MSG_INFO, stderr, "%s = %g\n", var->name, *((double *) var->variable));
            break;
        case JIBAL_CONFIG_VAR_UNIT:
            jabs_message(MSG_INFO, stderr, "%s = %g\n", var->name, *((double *) var->variable));
            break;
        case JIBAL_CONFIG_VAR_OPTION:
            jabs_message(MSG_INFO, stderr, "%s = %s\n", var->name,
                         jibal_option_get_string(var->option_list, *((int *) var->variable)));
            break;
        case JIBAL_CONFIG_VAR_SIZE:
            jabs_message(MSG_INFO, stderr, "%s = %zu\n", var->name, *((size_t *) var->variable));
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
        jabs_message(MSG_ERROR, stderr, "Variable %s is not boolean.\n", var->name);
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
        jabs_message(MSG_ERROR, stderr, "Variable %s is not boolean.\n", var->name);
    }
    return 0; /* Number of arguments */
}

script_command_status script_execute_command(script_session *s, const char *cmd) {
    int argc = 0;
    script_command_status status;
    if(!s) {
        jabs_message(MSG_ERROR, stderr, "Session has not been initialized.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    char **argv = string_to_argv(cmd, &argc);
    if(!argv) {
        jabs_message(MSG_ERROR, stderr, "Something went wrong in parsing arguments.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc) {
        status = script_execute_command_argv(s, s->commands, argc, argv); /* Note that s->file_depth may be altered (e.g. by script_load_script() */
    } else {
        status = SCRIPT_COMMAND_SUCCESS; /* Doing nothing successfully */
    }
    argv_free(argv, argc);
    return status;
}

script_command_status script_execute_command_argv(script_session *s, const script_command *commands, int argc, char **argv) {
    if(!s || !commands || !argv)
        return SCRIPT_COMMAND_FAILURE;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "No arguments given.\n");
        return SCRIPT_COMMAND_NOT_FOUND;
    }
    const script_command *cmds = commands;
    const script_command *c_parent = NULL;
    while(argc && cmds) { /* Arguments and subcommands remain. Try to find the right one, if possible. */
#ifdef DEBUG
        fprintf(stderr, "Top level script_execute_command_argv() loop, %i arguments remain (start with %s).\n", argc, argv[0]);
#endif
        const script_command *c = script_command_find(cmds, argv[0]);
        if(!c) {
#ifdef DEBUG
            fprintf(stderr, "Debug: Didn't find command %s.\n", argv[0]);
#endif
            script_command_not_found(argv[0], c_parent);
            return SCRIPT_COMMAND_NOT_FOUND;
        }
        while(c) { /* Subcommand found */
            argc--;
            argv++;
#ifdef DEBUG
            fprintf(stderr, "Debug: Found command %s.\n", c->name);
#endif
            if(c->f) {
#ifdef DEBUG
                fprintf(stderr, "Debug: There is a function %p in command %s. Calling it with %i arguments.\n", (void *) c->f, c->name, argc);
#endif
                script_command_status status = c->f(s, argc, argv);
                if(status > 0) { /* Positive numbers indicate number of arguments consumed */
                    argc -= status;
                    argv += status;
                }
#ifdef DEBUG
                fprintf(stderr, "Debug: Command run, returned %i (%s). Number of arguments remaining: %i\n", status,
                        script_command_status_to_string(status), argc);
#endif
                if(status < 0 && status != SCRIPT_COMMAND_NOT_FOUND) { /* Command not found is acceptable, we try to find subcommands later. All other errors cause an immediate return. */
                    return status;
                }
                if(status == SCRIPT_COMMAND_NOT_FOUND && argc == 0) { /* Command not found, no arguments remain, show an error. */
                    script_command_not_found(NULL, c);
                    return status;
                }
            } else if(c->var) {
#ifdef DEBUG
                fprintf(stderr, "Debug: %s is a var.\n", c->name);
#endif
                if(!c_parent) {
                    jabs_message(MSG_ERROR, stderr,
                                 "Command/option \"%s\" is a variable, but there is no parent command at all. This is highly unusual.\n",
                                 c->name);
                    return SCRIPT_COMMAND_FAILURE;
                }
                if(!c_parent->f_var) {
                    jabs_message(MSG_ERROR, stderr,
                                 "Command/option \"%s\" is a variable, but there is no function to handle variables in parent command (\"%s\"). This is highly unusual.\n",
                                 c->name, c_parent->name);
                    return SCRIPT_COMMAND_FAILURE;
                }
#ifdef DEBUG
                fprintf(stderr, "Debug: Using function f_var()=%p in %s. %i arguments. Arguments start with: %s\n", (void *)c_parent->f_var, c_parent->name, argc, argc?argv[0]:"(no arguments)");
#endif
                script_command_status status = c_parent->f_var(s, c->var, argc, argv);
                if(status >= 0) { /* Positive numbers indicate number of arguments consumed */
                    argc -= status;
                    argv += status;
                } else {
                    return status;
                }
            } else if(c->val) {
                if(!c_parent) {
                    jabs_message(MSG_ERROR, stderr,
                                 "Command/option \"%s\" has a non-zero value, but there is no parent command at all. This is highly unusual.\n",
                                 c->name);
                    return SCRIPT_COMMAND_FAILURE;
                }
                if(!c_parent->f_val) {
                    jabs_message(MSG_ERROR, stderr,
                                 "Command/option \"%s\" has a non-zero value, but there is no function to handle values in parent command (\"%s\"). This is highly unusual.\n",
                                 c->name, c_parent->name);
                    return SCRIPT_COMMAND_FAILURE;
                }
#ifdef DEBUG
                fprintf(stderr, "Debug: Using function f_val()=%p in %s. %i arguments. Arguments start with: %s\n", (void *)c_parent->f_val, c_parent->name, argc, argc?argv[0]:"(no arguments)");
#endif
                script_command_status status = c_parent->f_val(s, c->val, argc, argv);
                if(status >= 0) { /* Positive numbers indicate number of arguments consumed */
                    argc -= status;
                    argv += status;
                } else {
                    return status;
                }
            } else if(!argc) {
                script_command_not_found(NULL, c); /* No function, no nothing, no arguments. */
                return SCRIPT_COMMAND_NOT_FOUND;
            }
            if(c->subcommands) {
                cmds = c->subcommands;
                c_parent = c;
                break;
            } else {
#ifdef DEBUG
                fprintf(stderr,"Debug: there are no subcommands in \"%s\" (this is not an error). There is a val as always: %i.\n", c->name, c->val);
#endif
                c = NULL; /* Moving back to upper level. */
            }
        }
    }
    if(argc) {
        jabs_message(MSG_ERROR, stderr, "Debug: Didn't find command or an error with \"%s\" (%i arguments remain).\n", argv[0], argc);
        return SCRIPT_COMMAND_NOT_FOUND;
    }
    return SCRIPT_COMMAND_SUCCESS;
}

int command_compare(const void *a, const void *b) {
    const script_command *o_a = (const script_command *) a;
    const script_command *o_b = (const script_command *) b;
    return strcmp(o_a->name, o_b->name);
}

char *script_commands_list_matches(const script_command *commands, const char *str) {
    size_t len = strlen(str);
    size_t n = 0;
    int found = 0;
    for(const script_command *c = commands; c; c = c->next) {
        if(strncmp(c->name, str, len) == 0) { /* At least partial match */
            n += strlen(c->name) + 1;
            found++;
        }
    }
    n++;
    char *out = malloc(sizeof(char) * n);
    out[0] = '\0';
    for(const script_command *c = commands; c; c = c->next) {
        if(strncmp(c->name, str, len) == 0) { /* At least partial match */
            strcat(out, c->name);
            found--;
            if(found) {
                strcat(out, " ");
            }
        }
    }
    return out;
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

script_command *script_command_new(const char *name, const char *help_text, int val, script_command_status (*f)(struct script_session *, int, char *const *)) {
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

int script_command_set_var(script_command *c, jibal_config_var_type type, void *variable, const jibal_option *option_list) {
    if(!c->var) {
        c->var = malloc(sizeof(jibal_config_var));
    }
    c->var->variable = variable;
    c->var->type = type;
    c->var->name = c->name; /* The pointer is shared, so "var" doesn't get its own */
    c->var->option_list = option_list;
    c->f = NULL; /* These guys can't coexist */
    return EXIT_SUCCESS;
}

void script_command_free(script_command *c) {
    if(!c)
        return;
#ifdef DEBUG_VERBOSE
    fprintf(stderr, "Freeing command \"%s\" (%p)\n", c->name, (void *)c);
#endif
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
#ifdef DEBUG
    fprintf(stderr, "Sorting command list (head %p = %s)\n", (void *)head, head->name);
#endif
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
#ifdef WARN_ON_DUMMY_COMMANDS
    if(!c_new->subcommands && !c_new->f && !c_new->var && !c_new->val) { /* Everything needs to be set before calling this function to avoid this warning. */
        jabs_message(MSG_WARNING, stderr, "Warning: \"%s\" commands/option does not have subcommands and doesn't define a function, variable or a value!\n", c_new->name);
    }
#endif
    if(*head == NULL) {
        *head = c_new;
        return;
    }
    if(c_new->val) { /* Check if vals are multiply defined, because that would suck. */
        for(const script_command *c = *head; c; c = c->next) {
            if(c->val == c_new->val) {
                jabs_message(MSG_WARNING, stderr, "Warning: value %i (='%c') defined both in \"%s\" and \"%s\" commands/options!\n", c->val, c->val, c->name, c_new->name);
            }
        }
    }
    script_command *tail = script_command_list_find_tail(*head);
    if(tail) {
        tail->next = c_new;
    }
}

script_command *script_command_list_from_command_array(const script_command *commands) {
    if(!commands)
        return NULL;
    script_command *head = NULL;
    script_command *tail = NULL;
    for(const script_command *c = commands; c->name != NULL; c++) {
        script_command *c_new = script_command_new(c->name, c->help_text, c->val, c->f);
        if(c->var) {
            script_command_set_var(c_new, c->var->type, c->var->variable, c->var->option_list);
        }
        if(!head) {
            head = c_new;
            tail = c_new;
            continue;
        }
        tail->next = c_new;
        tail = c_new;
    }
    return head;
}

script_command *script_command_list_from_vars_array(const jibal_config_var *vars, jibal_config_var_type type) {
    script_command *head = NULL;
    for(const jibal_config_var *var = vars; var->type != 0; var++) {
        if(type != 0 && var->type != type) { /* Restrict by type */
            continue;
        }
        script_command *c = script_command_new(var->name, jibal_config_var_type_name(var->type), 0, NULL);
        script_command_set_var(c, var->type, var->variable, var->option_list);
        script_command_list_add_command(&head, c);
    }
    return head;
}

script_command *script_commands_create(struct script_session *s) {
    fit_data *fit = s->fit;
    simulation *sim = fit->sim;
    script_command *head = NULL;


    script_command *c;
    script_command *c_help = script_command_new("help", "Help.", 0, &script_help);
    script_command_list_add_command(&head, c_help);
    script_command_list_add_command(&c_help->subcommands, script_command_new("commands", "List of commands.", 0, &script_help_commands));
    script_command_list_add_command(&c_help->subcommands, script_command_new("version", "Help on (show) version.", 0, &script_help_version));

#ifdef JABS_PLUGINS
    script_command *c_identify = script_command_new("identify", "Identify something.", 0, NULL);
    script_command_list_add_command(&head, c_identify);
    script_command_list_add_command(&c_identify->subcommands, script_command_new("plugin", "Identify plugin.", 0, &script_identify_plugin));
#endif

    script_command *c_set = script_command_new("set", "Set something.", 0, NULL);
    c_set->f_var = &script_set_var;
    script_command_list_add_command(&head, c_set);
    script_command_list_add_command(&c_set->subcommands, script_command_new("aperture", "Set aperture.", 0, &script_set_aperture));

    script_command *c_detector = script_command_new("detector", "Set detector properties.", 0, &script_set_detector);
    c_detector->f_val = &script_set_detector_val;
    script_command_list_add_command(&c_set->subcommands, c_detector);
    script_command_list_add_command(&c_detector->subcommands, script_command_new("aperture", "Set detector aperture.", 0, &script_set_detector_aperture));
    script_command *c_calibration = script_command_new("calibration", "Set calibration.", 0, &script_set_detector_calibration);
    c_calibration->f_val = &script_set_detector_calibration_val;
    script_command_list_add_command(&c_detector->subcommands, c_calibration);
    script_command_list_add_command(&c_calibration->subcommands, script_command_new("linear", "Set the calibration to be linear (default).", 'L', NULL));
    script_command_list_add_command(&c_calibration->subcommands, script_command_new("slope", "Set the slope of a linear calibration.", 's', NULL));
    script_command_list_add_command(&c_calibration->subcommands, script_command_new("offset", "Set the offset of a linear calibration.", 'o', NULL));
    script_command_list_add_command(&c_calibration->subcommands, script_command_new("resolution", "Set the resolution.", 'r', NULL));
    script_command_list_add_command(&c_calibration->subcommands, script_command_new("poly", "Set the calibration to be a polynomial.", 0, &script_set_detector_calibration_poly));

    script_command_list_add_command(&c_detector->subcommands, script_command_new("column", "Set column number (for data input).", 'c', NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("channels", "Set number of channels.", 'h', NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("compress", "Set compress (summing of channels).", 'C', NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("distance", "Set detector distance from target.", 'd', NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("foil", "Set detector foil.", 0, &script_set_detector_foil));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("type", "Set detector type.", 't', NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("slope", "Set detector calibration slope.", 'S', NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("offset", "Set detector calibration offset.", 'O', NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("resolution", "Set detector resolution (FWHM).", 'r', NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("solid", "Set detector solid angle.", 's', NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("theta", "Set detector (scattering) angle.", 'T', NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("length", "Set detector length (for ToF).", 'l', NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("beta", "Set exit angle (angle of ion in sample) for this detector.", 'b', NULL));
    script_command_list_add_command(&c_detector->subcommands, script_command_new("phi", "Set detector azimuth angle, 0 = IBM, 90 deg = Cornell.", 'p', NULL));

    script_command_list_add_command(&c_set->subcommands, script_command_new("ion", "Set incident ion (isotope).", 0, &script_set_ion));
    script_command *c_set_fit = script_command_new("fit", "Set fit related things.", 0, NULL);
    c_set_fit->f_val = &script_set_fit_val;
    script_command_list_add_command(&c_set_fit->subcommands, script_command_new("normal", "Normal two-phase fitting.", 'n', NULL));
    script_command_list_add_command(&c_set_fit->subcommands, script_command_new("slow", "One phase fitting (slow phase only).", 's', NULL));
    script_command_list_add_command(&c_set_fit->subcommands, script_command_new("fast", "One phase fitting (fast phase only).", 'f', NULL));
    script_command_list_add_command(&c_set->subcommands, c_set_fit);

    script_command_list_add_command(&c_set->subcommands, script_command_new("sample", "Set sample.", 0, &script_set_sample));

    script_command *c_set_simulation = script_command_new("simulation", "Set simulation related things.", 0, NULL);
    c_set_simulation->f_val = &script_set_simulation_val;
    script_command_list_add_command(&c_set_simulation->subcommands, script_command_new("defaults", "Set default calculation parameters.", 'd', NULL));
    script_command_list_add_command(&c_set_simulation->subcommands, script_command_new("brisk", "Set slightly faster calculation parameters.", 'b', NULL));
    script_command_list_add_command(&c_set_simulation->subcommands, script_command_new("fast", "Set fast calculation parameters.", 'f', NULL));
    script_command_list_add_command(&c_set_simulation->subcommands, script_command_new("accurate", "Set the most accurate calculation parameters.", 'a', NULL));
    script_command_list_add_command(&c_set_simulation->subcommands, script_command_new("improved", "Set more accurate calculation parameters.", 'i', NULL));
    script_command_list_add_command(&c_set->subcommands, c_set_simulation);

    script_command_list_add_command(&c_set->subcommands, script_command_new("stopping", "Set (assign) stopping or straggling.", 0, &script_set_stopping));
    script_command_list_add_command(&c_set->subcommands, script_command_new("variable", "Set a variable.", 0, NULL));

    const jibal_config_var vars[] = {
            {JIBAL_CONFIG_VAR_UNIT,   "fluence",                     &sim->fluence,                             NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "energy",                      &sim->beam_E,                              NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "energy_broad",                &sim->beam_E_broad,                        NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "emin",                        &sim->emin,                                NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "alpha",                       &sim->sample_theta,                        NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "phi",                         &sim->sample_phi,                          NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "channeling",                  &sim->channeling_offset,                   NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "channeling_slope",            &sim->channeling_slope,                    NULL},
            {JIBAL_CONFIG_VAR_STRING, "output",                      &s->output_filename,                       NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "erd",                         &sim->erd,                                 NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "rbs",                         &sim->rbs,                                 NULL},
            {JIBAL_CONFIG_VAR_SIZE,   "maxiter",                     &fit->n_iters_max,                         NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "xtolerance",                  &fit->xtol,                                NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "gtolerance",                  &fit->gtol,                                NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "ftolerance",                  &fit->ftol,                                NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "chisq_tolerance",             &fit->chisq_tol,                           NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "chisq_fast_tolerance",        &fit->chisq_fast_tol,                      NULL},
            {JIBAL_CONFIG_VAR_SIZE,   "depthsteps_max",              &sim->params->depthsteps_max,              NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "ds",                          &sim->params->ds,                          NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "rk4",                         &sim->params->rk4,                         NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "sigmas_cutoff",               &sim->params->sigmas_cutoff,               NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "stop_step_incident",          &sim->params->stop_step_incident,          NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "stop_step_exiting",           &sim->params->stop_step_exiting,           NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "stop_step_fudge",             &sim->params->stop_step_fudge_factor,      NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "stop_step_min",               &sim->params->stop_step_min,               NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "stop_step_add",               &sim->params->stop_step_add,               NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "nucl_stop_accurate",          &sim->params->nucl_stop_accurate,          NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "mean_conc_and_energy",        &sim->params->mean_conc_and_energy,        NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "geostragg",                   &sim->params->geostragg,                   NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "beta_manual",                 &sim->params->beta_manual,                 NULL},
            {JIBAL_CONFIG_VAR_SIZE,   "cs_n_steps",                  &sim->params->cs_n_steps,                  NULL},
            {JIBAL_CONFIG_VAR_SIZE,   "cs_n_stragg_steps",           &sim->params->cs_n_stragg_steps,           NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "gaussian_accurate",           &sim->params->gaussian_accurate,           NULL},
            {JIBAL_CONFIG_VAR_SIZE,   "int_cs_max_intervals",        &sim->params->int_cs_max_intervals,        NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "int_cs_accuracy",             &sim->params->int_cs_accuracy,             NULL},
            {JIBAL_CONFIG_VAR_SIZE,   "int_cs_stragg_max_intervals", &sim->params->int_cs_stragg_max_intervals, NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "int_cs_stragg_accuracy",      &sim->params->int_cs_stragg_accuracy,      NULL},
            {JIBAL_CONFIG_VAR_NONE, NULL, NULL,                                                                 NULL}
    };
    c = script_command_list_from_vars_array(vars, 0);
    script_command_list_add_command(&c_set->subcommands, c);

    script_command *c_load = script_command_new("load", "Load something.", 0, NULL);
    script_command_list_add_command(&head, c_load);
    script_command_list_add_command(&c_load->subcommands, script_command_new("experimental", "Load an experimental spectrum.", 0, &script_load_experimental));
    script_command_list_add_command(&c_load->subcommands, script_command_new("script", "Load (run) a script.", 0, &script_load_script));
    script_command_list_add_command(&c_load->subcommands, script_command_new("sample", "Load a sample.", 0, &script_load_sample));
    script_command *c_reaction = script_command_new("reaction", "Load a reaction from R33 file.", 0, &script_load_reaction);
    script_command_list_add_command(&c_load->subcommands, c_reaction);
#ifdef JABS_PLUGINS
    script_command_list_add_command(&c_reaction->subcommands, script_command_new("plugin", "Load a reaction from a plugin.", 0, &script_load_reaction_plugin));
#endif
    script_command_list_add_command(&c_load->subcommands, script_command_new("roughness", "Load layer thickness table (roughness) from a file.", 0, &script_load_roughness));

    script_command *c_show = script_command_new("show", "Show information on things.", 0, NULL);
    script_command_list_add_command(&head, c_show);
    script_command_list_add_command(&c_show->subcommands, script_command_new("aperture", "Show aperture.", 0, &script_show_aperture));
    script_command_list_add_command(&c_show->subcommands, script_command_new("detector", "Show detector.", 0, &script_show_detector));

    script_command *c_fit = script_command_new("fit", "Show fit results.", 0, &script_show_fit);
    script_command_list_add_command(&c_show->subcommands, c_fit);
    script_command_list_add_command(&c_fit->subcommands, script_command_new("variables", "Show possible fit variables.", 0, &script_show_fit_variables));
    script_command_list_add_command(&c_fit->subcommands, script_command_new("ranges", "Show fit ranges.", 0, &script_show_fit_ranges));

    script_command_list_add_command(&c_show->subcommands, script_command_new("reactions", "Show reactions.", 0, &script_show_reactions));
    script_command_list_add_command(&c_show->subcommands, script_command_new("sample", "Show sample.", 0, &script_show_sample));
    script_command_list_add_command(&c_show->subcommands, script_command_new("simulation", "Show simulation.", 0, &script_show_simulation));
    script_command_list_add_command(&c_show->subcommands, script_command_new("stopping", "Show stopping (GSTO) assignments.", 0, &script_show_stopping));

    script_command *c_show_variable = script_command_new("variable", "Show variable.", 0, NULL);
    c_show_variable->f_var = script_show_var;
    script_command_list_add_command(&c_show->subcommands, c_show_variable);
    c = script_command_list_from_vars_array(vars, 0);
    script_command_list_add_command(&c_show_variable->subcommands, c);

    script_command_list_add_command(&head, script_command_new("exit", "Exit.", 0, &script_exit));

    c = script_command_new("save", "Save something.", 0, NULL);
    script_command_list_add_command(&head, c);
    script_command_list_add_command(&c->subcommands, script_command_new("bricks", "Save bricks.", 0, &script_save_bricks));
    script_command_list_add_command(&c->subcommands, script_command_new("calibrations", "Save detector calibrations.", 0, &script_save_calibrations));
    script_command_list_add_command(&c->subcommands, script_command_new("sample", "Save sample.", 0, &script_save_sample));
    script_command_list_add_command(&c->subcommands, script_command_new("spectra", "Save spectra.", 0, &script_save_spectra));

    c = script_command_new("test", "Test something.", 0, NULL);
    script_command_list_add_command(&head, c);
    script_command_list_add_command(&c->subcommands, script_command_new("file", "Test simulated spectrum against a reference.", 0, &script_test_file));
    script_command_list_add_command(&c->subcommands, script_command_new("roi", "Test ROI.", 0, &script_test_roi));

    c = script_command_new("add", "Add something.", 0, NULL);
    script_command_list_add_command(&head, c);
    script_command_list_add_command(&c->subcommands, script_command_new("detector", "Add a detector.", 0, &script_add_detector));

    script_command *c_add_fit = script_command_new("fit", "Add something related to fit.", 0, NULL);
    script_command_list_add_command(&c->subcommands, c_add_fit);
    script_command_list_add_command(&c_add_fit->subcommands, script_command_new("range", "Add a fit range", 0, &script_add_fit_range));

    script_command_list_add_command(&c->subcommands, script_command_new("reaction", "Add a reaction.", 0, &script_add_reaction));
    script_command_list_add_command(&c->subcommands, script_command_new("reactions", "Add reactions (of some type).", 0, &script_add_reactions));

    c = script_command_new("remove", "Remove something.", 0, NULL);
    script_command_list_add_command(&head, c);
    script_command_list_add_command(&c->subcommands, script_command_new("reaction", "Remove reaction.", 0, &script_remove_reaction));

    c = script_command_new("reset", "Reset something (or everything).", 0, &script_reset);
    script_command_list_add_command(&head, c);
    script_command_list_add_command(&c->subcommands, script_command_new("detectors", "Reset detectors.", 0, &script_reset_detectors));
    script_command_list_add_command(&c->subcommands, script_command_new("experimental", "Reset experimental spectra.", 0, &script_reset_experimental));
    script_command_list_add_command(&c->subcommands, script_command_new("fit", "Reset fit (ranges).", 0, &script_reset_fit));
    script_command_list_add_command(&c->subcommands, script_command_new("reactions", "Reset reactions.", 0, &script_reset_reactions));
    script_command_list_add_command(&c->subcommands, script_command_new("sample", "Reset sample.", 0, &script_reset_sample));
    script_command_list_add_command(&c->subcommands, script_command_new("stopping", "Reset stopping assignments.", 0, &script_reset_stopping));

    c = script_command_new("fit", "Do a fit.", 0, script_fit);
    script_command_list_add_command(&head, c);

    script_command *c_enable = script_command_new("enable", "Set boolean variable to true.", 0, NULL);
    c_enable->f_var = &script_enable_var;
    script_command_list_add_command(&head, c_enable);
    c = script_command_list_from_vars_array(vars, JIBAL_CONFIG_VAR_BOOL);
    script_command_list_add_command(&c_enable->subcommands, c);

    script_command *c_disable = script_command_new("disable", "Set boolean variable to true.", 0, NULL);
    c_disable->f_var = &script_disable_var;
    script_command_list_add_command(&head, c_disable);
    c = script_command_list_from_vars_array(vars, JIBAL_CONFIG_VAR_BOOL);
    script_command_list_add_command(&c_disable->subcommands, c);

    c = script_command_new("roi", "Show information from a region of interest.", 0, script_roi);
    script_command_list_add_command(&head, c);

    c = script_command_new("simulate", "Run a simulation.", 0, script_simulate);
    script_command_list_add_command(&head, c);

    script_command *c_split = script_command_new("split", "Split something.", 0, NULL);
    script_command_list_add_command(&head, c_split);
    script_command *c_split_sample = script_command_new("sample", "Split something sample related.", 0, NULL);
    script_command_list_add_command(&c_split->subcommands, c_split_sample);
    script_command_list_add_command(&c_split_sample->subcommands, script_command_new("elements", "Split materials down to their constituent elements.", 0, script_split_sample_elements));

    script_command_list_add_command(&head, script_command_new("cwd", "Display current working directory.", 0, script_cwd));
    script_command_list_add_command(&head, script_command_new("pwd", "Display current working directory.", 0, script_cwd));
    script_command_list_add_command(&head, script_command_new("cd", "Change current working directory.", 0, script_cd));
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
#ifdef DEBUG
            fprintf(stderr, "Sorting all commands, level %zu, subcommands of %s\n", i, c->name);
#endif
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

void script_commands_print(FILE *f, const struct script_command *commands) {
    if(!commands)
        return;
    for(const struct script_command *c = commands; c; c = c->next) {
        if(!c->help_text)
            continue;
        jabs_message(MSG_INFO, f, " %20s    %s\n", c->name, c->help_text);
    }
}

size_t script_commands_size(const script_command *commands) {
    if(!commands)
        return 0;
    size_t n = 0;
    for(const struct script_command *c = commands; c; c++) {
        n++;
    }
#ifdef DEBUG
    fprintf(stderr, "Commands size is %zu (%p).\n", n, (void *)commands);
#endif
    return n;
}

void script_print_command_tree(FILE *f, const struct script_command *commands) {
    const struct script_command *stack[SCRIPT_COMMANDS_NESTED_MAX];
    const struct script_command *c;
    stack[0] = commands;
    size_t i = 0;
    c = stack[0];
    while(c) {
#ifdef DEBUG
        if(TRUE) {
#else
        if(c->f || c->var || c->val) { /* If none of these is set, we shouldn't print the command name at all */
#endif
            for(size_t j = 0; j < i; j++) {
                jabs_message(MSG_INFO, f, "%s ", stack[j]->name);
            }
            jabs_message(MSG_INFO, f, "%s\n", c->name);
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
        jabs_message(MSG_ERROR, stderr, "Usage: load script [file]\n");
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
        jabs_message(MSG_ERROR, stderr, "Usage: load sample [file]\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    sample_model *sm = sample_model_from_file(fit->jibal, argv[0]);
    if(!sm) {
        jabs_message(MSG_INFO, stderr, "Sample load from \"%s\" failed.\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    sample_model_free(fit->sm);
    fit->sm = sm;
    if(s->fit->sim->n_reactions > 0) {
        jabs_message(MSG_WARNING, stderr, "Reactions were reset automatically, since the sample was changed.\n");
        sim_reactions_free(fit->sim);
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
        jabs_message(MSG_ERROR, stderr,
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
                const gsl_histogram *h = fit->exp[i_det];
                jabs_message(MSG_VERBOSE, stderr, "Detector %zu: experimental spectrum with %zu channels loaded.\n", i_det + 1, h?h->n:0);
            }
        }
    }
    argc--;
    argv++;
    return argc_orig - argc; /* Number of arguments */
}

script_command_status script_load_reaction(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: load reaction <file>\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(strcmp(argv[0], "plugin") == 0) {
        return 0;
    }
    if(sim_reactions_add_r33(fit->sim, fit->jibal->isotopes, argv[0])) {
        return SCRIPT_COMMAND_FAILURE;
    } else {
        return 1;
    }
}

script_command_status script_load_roughness(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    sample_model *sm = fit->sm;
    if(argc < 2) {
        jabs_message(MSG_ERROR, stderr, "Usage: load roughness <layer> <file>\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(!sm) {
        jabs_message(MSG_ERROR, stderr, "No sample has been set, can't use this command yet.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    size_t n = strtoull(argv[0], NULL, 10);
    if(n == 0 || n > sm->n_ranges) {
        jabs_message(MSG_ERROR, stderr, "Layer number must be between 1 and %zu (number of ranges in current sample).\n", sm->n_ranges);
        return SCRIPT_COMMAND_FAILURE;
    }
    sample_range *range = &(sm->ranges[n - 1]);
    if(roughness_set_from_file(&range->rough, argv[1])) {
        jabs_message(MSG_ERROR, stderr, "Setting roughness from file failed.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    range->x = thickness_probability_table_areal_density(range->rough.file->tpd); /* Update thickness of layer to correspond the average areal density */
    if(range->rough.file && range->rough.file->tpd) {
        jabs_message(MSG_INFO, stderr, "Layer %zu: roughness from file \"%s\", containing %zu data points, loaded.\n", n + 1, range->rough.file->filename, range->rough.file->tpd->n);
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
    fit_data_fit_ranges_free(s->fit);
    fit_params_free(s->fit->fit_params);
    s->fit->fit_params = NULL;
    fit_data_workspaces_free(s->fit);
    fit_data_histo_sum_free(s->fit);
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
    struct fit_data *fit = s->fit;
    fit_data_exp_free(s->fit);
    fit->exp = calloc(fit->sim->n_det, sizeof(gsl_histogram *));
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
#ifdef DEBUG
    fprintf(stderr, "Resetting everything!\n");
#endif
    fit_data_fit_ranges_free(fit);
    fit_params_free(fit->fit_params);
    fit->fit_params = NULL;
    fit_data_exp_free(fit);
    fit_data_workspaces_free(fit);
    sample_model_free(fit->sm);
    fit->sm = NULL;
    sim_free(fit->sim);
    fit_data_free(s->fit);
    s->fit = fit_data_new(s->jibal, sim_init(s->jibal));
    if(!s->fit) {
        jabs_message(MSG_ERROR, stderr, "Reset fails due to an allocation issue. This should not happen.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    script_commands_free(s->commands);
    s->commands = script_commands_create(s);
    jibal_gsto_assign_clear_all(s->jibal->gsto);
    return 0;
}

script_command_status script_show_sample(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    struct fit_data *fit = s->fit;
    if(!fit->sm) {
        jabs_message(MSG_WARNING, stderr, "No sample has been set.\n");
        return 0;
    }
#ifdef DEBUG
    char *sample_str = sample_model_to_string(fit->sm);
        fprintf(stderr, "Sample: %s\n", sample_str);
        free(sample_str);
#endif
    jabs_message(MSG_INFO, stderr, "Sample model (use \"save sample\" to save):\n");
    if(sample_model_print(NULL, fit->sm)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    if(fit->sm->type == SAMPLE_MODEL_LAYERED) {
        jabs_message(MSG_INFO, stderr, "\nLayers:\n");
        sample *sample = sample_from_sample_model(fit->sm);
        sample_print_thicknesses(NULL, sample);
        sample_free(sample);
    }
    return 0;
}

script_command_status script_show_simulation(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    sim_print(s->fit->sim);
    return 0;
}

script_command_status script_show_stopping(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    jibal_gsto *gsto = s->jibal->gsto;
    int n = 0;
    jabs_message(MSG_INFO, stderr, "List of assigned stopping and straggling files:\n");
    for(int Z1 = 1; Z1 <= gsto->Z1_max; Z1++) {
        for(int Z2 = 1; Z2 <= gsto->Z2_max; Z2++) {
            gsto_file_t *file_sto = jibal_gsto_get_assigned_file(gsto, GSTO_STO_ELE, Z1, Z2);
            gsto_file_t *file_stg = jibal_gsto_get_assigned_file(gsto, GSTO_STO_STRAGG, Z1, Z2);
            if(!file_sto && !file_stg) {
                continue; /* Nothing to do, nothing assigned */
            }
            jabs_message(MSG_INFO, stderr, "  Z1=%i (%s), Z2=%i (%s): ", Z1, jibal_element_name(gsto->elements, Z1), Z2,
                         jibal_element_name(gsto->elements, Z2));
            if(file_sto) {
                jabs_message(MSG_INFO, stderr, "Stopping file %s.", file_sto->name);
                n++;
            }
            if(file_stg) {
                jabs_message(MSG_INFO, stderr, "%sStraggling file %s.", file_sto ? " " : "", file_stg->name);
                n++;
            }
            jabs_message(MSG_INFO, stderr, "\n");
        }
    }
    jabs_message(MSG_INFO, stderr, "\nTotal of %i assignments.\n", n);
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_show_fit(script_session *s, int argc, char *const *argv) {
    (void) argv;
    if(argc == 0) {
        if(!s->fit->fit_params || s->fit->fit_params->n_active == 0) {
            jabs_message(MSG_INFO, stderr, "No fit results to show.\n");
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
    fit_params_print(p_all, FALSE, pattern);
    fit_params_free(p_all);
    return argc_orig - argc;
}

script_command_status script_show_fit_ranges(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    fit_data_print(stderr, s->fit);
    return 0;
}

script_command_status script_show_aperture(struct script_session *s, int argc, char *const *argv) {
    (void) argv;
    struct fit_data *fit = s->fit;
    const int argc_orig = argc;
    char *aperture_str = aperture_to_string(fit->sim->beam_aperture);
    jabs_message(MSG_INFO, stderr, "aperture %s\n", aperture_str);
    free(aperture_str);
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
        jabs_message(MSG_INFO, stderr, "No detectors have been defined.\n");
    } else if(argc != argc_orig || fit->sim->n_det == 1) { /* Show information on particular detector if a number was given or if we only have one detector. */
        if(detector_print(s->jibal, NULL, sim_det(fit->sim, i_det))) {
            jabs_message(MSG_ERROR, stderr, "No detectors set or other error.\n");
        }
    } else {
        jabs_message(MSG_INFO, stderr, "  # | col | theta |  phi  |  solid  |   type   | calibration\n");
        for(size_t i = 0; i < fit->sim->n_det; i++) {
            detector *det = sim_det(fit->sim, i);
            if(!det)
                continue;
            char *calib_str = calibration_to_string(det->calibration);
            jabs_message(MSG_INFO, stderr, "%3zu | %3zu | %5.1lf | %5.1lf | %7.3lf | %8s | %s\n",
                         i + 1, det->column, det->theta / C_DEG, det->phi / C_DEG, det->solid / C_MSR, detector_type_name(det), calib_str);
            free(calib_str);
        }
        jabs_message(MSG_INFO, stderr, "Use 'show detector <number>' to get more information on a particular detector.\n");
    }
    return argc_orig - argc; /* Number of arguments */
}

script_command_status script_show_reactions(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    struct fit_data *fit = s->fit;
    if(fit->sim->n_reactions == 0) {
        jabs_message(MSG_INFO, stderr, "No reactions.\n");
    }
    reactions_print(stderr, fit->sim->reactions, fit->sim->n_reactions);
    return 0;
}

script_command_status script_set_ion(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: set ion <isotope>\nExample: set ion 4He\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    const jibal_isotope *isotope = jibal_isotope_find(fit->jibal->isotopes, argv[0], 0, 0);
    if(!isotope) {
        jabs_message(MSG_ERROR, stderr, "No such isotope: %s\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    s->fit->sim->beam_isotope = isotope;
    return 1;
}

script_command_status script_set_aperture(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    const int argc_orig = argc;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: set aperture <type> {width|height|diameter} ...\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    fit->sim->beam_aperture = aperture_set_from_argv(s->jibal, fit->sim->beam_aperture, &argc, &argv);
    if(!fit->sim->beam_aperture) {
        jabs_message(MSG_ERROR, stderr, "Aperture could not be parsed.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc) {
        jabs_message(MSG_ERROR, stderr, "Unexpected extra arguments (%i), starting with %s.\n", argc, argv[0]);
        if(fit->sim->beam_aperture->type == APERTURE_NONE) {
            jabs_message(MSG_INFO, stderr, "Aperture type not defined. Allowed types:");
            for(const jibal_option *o = aperture_option; o->s; o++) {
                jabs_message(MSG_INFO, stderr, " %s", o->s);
            }
            jabs_message(MSG_INFO, stderr, "\n");
        } else {
            jabs_message(MSG_ERROR, stderr, "Aperture type was %s. Aperture removed.\n", aperture_name(fit->sim->beam_aperture));
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
        jabs_message(MSG_ERROR, stderr, "Usage: set detector aperture <type> {width|height|diameter} ...\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    detector *det = sim_det(s->fit->sim, s->i_det_active);
    if(!det) {
        jabs_message(MSG_ERROR, stderr, "No detector(s)\n");
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
#ifdef DEBUG
    fprintf(stderr, "Active Z is now %i\n", s->Z_active);
#endif
    return argc_orig - argc;
}

script_command_status script_set_detector_foil(struct script_session *s, int argc, char *const *argv) {
    const int argc_orig = argc;
    detector *det = sim_det(s->fit->sim, s->i_det_active);
    if(!det) {
        jabs_message(MSG_ERROR, stderr, "No detector(s)\n");
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
        jabs_message(MSG_ERROR, stderr, "No detector(s)\n");
        return SCRIPT_COMMAND_SUCCESS;
    }
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr,
                     "Usage: set detector calibration poly <n> <p_0> <p_1> ... <p_(n+1)>\nExample: set calibration poly 2 10keV 1.0keV 0.001keV\nThe example sets a second degree (quadratic 3 parameters) polynomial calibration.\nThe first parameter (p_0) is the constant term.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    size_t n = strtoull(argv[0], NULL, 10);
    argc--;
    argv++;
    if(argc < (int) (n + 1)) {
        jabs_message(MSG_ERROR, stderr, "Not enough parameters for a %zu degree polynomial. Expected %zu.\n", n, n + 1);
        return SCRIPT_COMMAND_FAILURE;
    }
    calibration *c = calibration_init_poly(n);
    calibration_copy_params(c, detector_get_calibration(det, s->Z_active)); /* This, de facto, only copies resolution from old calibration, since the rest are overwritten very soon. */
    for(int i = 0; i <= (int) n; i++) {
        calibration_set_param(c, i, jibal_get_val(s->jibal->units, UNIT_TYPE_ANY, argv[0]));
        argc--;
        argv++;
    }
    detector_set_calibration_Z(s->jibal->config, det, c, s->Z_active);
    return argc_orig - argc;
}

script_command_status script_set_sample(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 2) {
        jabs_message(MSG_ERROR, stderr, "Usage: set sample {nosimplify} <sample description>\nExample: set sample TiO2 1000tfu Si 10000tfu\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    const int argc_orig = argc;
    sample_model *sm_new = sample_model_from_argv(fit->jibal, &argc, &argv);
    if(!sm_new) {
        jabs_message(MSG_WARNING, stderr, "Setting sample fails.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    int argc_consumed = argc_orig - argc;
    sample_model_free(fit->sm);
    fit->sm = sm_new;
    if(s->fit->sim->n_reactions > 0) {
        jabs_message(MSG_WARNING, stderr, "Reactions were reset automatically, since the sample was changed.\n");
        sim_reactions_free(fit->sim);
    }
    return argc_consumed;
}

script_command_status script_set_stopping(struct script_session *s, int argc, char *const *argv) {
    const int argc_orig = argc;
    if(argc < 3) {
        jabs_message(MSG_ERROR, stderr, "Usage: set stopping <file> <element or Z1> <element or Z2> {<element or Z2>}\n");
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
        jabs_message(MSG_ERROR, stderr, "Error: \"%s\" is not a file recognized by GSTO.\n", argv[0]);
        if(s->jibal->gsto->n_files) {
            jabs_message(MSG_INFO, stderr, "List of files:");
            for(size_t i_file = 0; i_file < s->jibal->gsto->n_files; i_file++) {
                jabs_message(MSG_INFO, stderr, " %s", s->jibal->gsto->files[i_file].name);
            }
            jabs_message(MSG_INFO, stderr, "\n");
        } else {
            jabs_message(MSG_INFO, stderr, "There are no GSTO files, they should be listed in this file: \"%s\"\n",
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
            jabs_message(MSG_ERROR, stderr, "Could not assign stopping for Z1 = %i in Z2 = %i to file \"%s\".\n",
                         first->Z, second->Z, file->name);
            return SCRIPT_COMMAND_FAILURE;
        }
    }
    return argc_orig - argc;
}

script_command_status script_test_file(struct script_session *s, int argc, char *const *argv) {
    const struct fit_data *fit = s->fit;
    const int argc_orig = argc;
    size_t i_det = 0;
    if(script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc < 3) {
        jabs_message(MSG_ERROR, stderr, "Usage: test file <range> <filename> <tolerance>\nTests if simulation is within tolerance to a reference file.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    const gsl_histogram *sim = fit_data_sim(fit, i_det);
    if(!sim) {
        jabs_message(MSG_ERROR, stderr, "No simulation.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    const detector *det = sim_det(fit->sim, i_det);
    roi r = {.i_det = i_det};
    if(fit_set_roi_from_string(&r, argv[0])) {
        jabs_message(MSG_ERROR, stderr, "Could not parse range.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    gsl_histogram *h_ref = spectrum_read_detector(argv[1], det);
    if(!h_ref) {
        jabs_message(MSG_ERROR, stderr, "Could not read spectrum.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    double tolerance = strtod(argv[2], NULL);
    double error;
    argc -= 3;
    argv += 3;
    int return_value = argc_orig - argc;
    if(spectrum_compare(sim, h_ref, r.low, r.high, &error)) { /* Failed, test fails */
        jabs_message(MSG_ERROR, stderr, "Test failed. Is range valid?\n");
        return_value = SCRIPT_COMMAND_FAILURE;
    } else {
        jabs_message(MSG_INFO, stderr, "Test of simulated spectrum to reference from %zu to %zu. Error %e.\n", r.low, r.high, error);
        if(error > tolerance) {
            jabs_message(MSG_ERROR, stderr, "Test failed.\n");
            return_value = SCRIPT_COMMAND_FAILURE;
        } else {
            jabs_message(MSG_ERROR, stderr, "Test passed.\n");
        }
    }
    gsl_histogram_free(h_ref);
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
        jabs_message(MSG_ERROR, stderr, "Usage: test roi {<detector>} {exp} <range> <sum> <tolerance>\nTests if simulated (or experimental) ROI sum is within relative tolerance.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    roi r = {.i_det = i_det};
    if(fit_set_roi_from_string(&r, argv[0])) {
        jabs_message(MSG_ERROR, stderr, "Could not parse range.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    double sum_ref = strtod(argv[1], NULL);
    double tolerance = strtod(argv[2], NULL);
    gsl_histogram *h = exp ? fit_data_exp(fit, i_det) : fit_data_sim(fit, i_det);
    if(!h) {
        jabs_message(MSG_ERROR, stderr, "No histogram.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    double sum = spectrum_roi(h, r.low, r.high);
    argc -= 3;
    argv += 3;
    double rel_err = fabs(1.0 - sum / sum_ref);
    jabs_message(MSG_INFO, stderr, "Test of ROI from %zu to %zu. Sum %g. Relative error %e.\n", r.low, r.high, sum, rel_err);
    if(rel_err < tolerance) {
        jabs_message(MSG_INFO, stderr, "Test passed.\n");
    } else {
        jabs_message(MSG_ERROR, stderr, "Test failed.\n");
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
        jabs_message(MSG_ERROR, stderr, "Usage: add reaction TYPE isotope\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    const int argc_orig = argc;
    reaction *r = reaction_make_from_argv(fit->jibal, fit->sim->beam_isotope, &argc, &argv);
    int argc_consumed = argc_orig - argc;
    if(!r) {
        jabs_message(MSG_ERROR, stderr, "Could not make a reaction based on given description.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(r->cs == JIBAL_CS_NONE) {
        jibal_cross_section_type cs = sim_cs(fit->sim, r->type);
        r->cs = cs;
        jabs_message(MSG_VERBOSE, stderr,
                     "Reaction cross section not given or not valid, assuming default for %s: %s.\n",
                     reaction_name(r), jibal_cs_types[cs].s);
    }
    if(sim_reactions_add_reaction(fit->sim, r)) {
        return SCRIPT_COMMAND_FAILURE;
    } else {
        return argc_consumed;
    }
}

script_command_status script_add_reactions(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    if(!fit->sm) {
        jabs_message(MSG_ERROR, stderr,
                     "Cannot add reactions before sample has been set (I need to know which reactions to add!).\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc > 0) {
        if(strcmp(argv[0], "RBS") == 0) {
            sim_reactions_add_auto(fit->sim, fit->sm, REACTION_RBS, sim_cs(fit->sim, REACTION_RBS));
        } else if(strcmp(argv[0], "RBS-") == 0) {
            sim_reactions_add_auto(fit->sim, fit->sm, REACTION_RBS_ALT, sim_cs(fit->sim, REACTION_RBS_ALT));
        } else if(strcmp(argv[0], "ERD") == 0) {
            sim_reactions_add_auto(fit->sim, fit->sm, REACTION_ERD, sim_cs(fit->sim, REACTION_ERD));
        } else {
            jabs_message(MSG_ERROR, stderr, "What is %s anyways? Did you mean \"RBS\" or \"ERD\"?\n", argv[0]);
            return SCRIPT_COMMAND_FAILURE;
        }
        return 1;
    }
    if(fit->sim->rbs) {
        sim_reactions_add_auto(fit->sim, fit->sm, REACTION_RBS,
                               sim_cs(fit->sim, REACTION_RBS)); /* TODO: loop over all detectors and add reactions that are possible (one reaction for all detectors) */
        sim_reactions_add_auto(fit->sim, fit->sm, REACTION_RBS_ALT,
                               sim_cs(fit->sim, REACTION_RBS_ALT));
    }

    if(sim_do_we_need_erd(fit->sim)) {
        sim_reactions_add_auto(fit->sim, fit->sm, REACTION_ERD, sim_cs(fit->sim, REACTION_ERD));
    }
    return 0;
}

script_command_status script_add_detector(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: add detector {default}\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    script_set_detector(s, argc - 1, argv + 1);
    detector *det;
    if(strcmp(argv[0], "default") == 0) {
        det = detector_default(NULL);
    } else {
        jabs_message(MSG_WARNING, stderr, "Adding other types of detectors except default (add detector default) is currently not supported.\n");
        return 0; /* No arguments consumed, no error */
    }
    if(fit_data_add_det(fit, det)) {
        return SCRIPT_COMMAND_FAILURE;
    } else {
        return 1;
    }
}

script_command_status script_add_fit_range(script_session *s, int argc, char *const *argv) {
    struct fit_data *fit = s->fit;
    size_t i_det = 0;
    const int argc_orig = argc;
    script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det);
    int n_ranges = 0;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr,
                     "Usage: add fit_range {detector} <range> {<range> <range> ...}\nExample: add fit_range 1 [400:900] [980:1200]\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    while(argc > 0) {
        roi r = {.i_det = i_det};
        if(fit_set_roi_from_string(&r, argv[0])) {
            if(n_ranges == 0) {
                jabs_message(MSG_ERROR, stderr, "No ranges added! Failed on parsing \"%s\".\n", argv[0]);
                return SCRIPT_COMMAND_FAILURE;
            } else {
                break;
            }
        }
        if(fit_data_fit_range_add(fit, &r)) {
            jabs_message(MSG_ERROR, stderr, "Range not valid! Failed on parsing \"%s\".\n", argv[0]);
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
        jabs_message(MSG_INFO, stderr, "Type help [topic] for information on a particular topic or \"help help\" for help on help.\n\n");

        return SCRIPT_COMMAND_NOT_FOUND;
    }

    int found = 0;
    for(const struct help_topic *t = topics; t->name != NULL; t++) {
        if(strcmp(t->name, argv[0]) == 0) {
            found++;
            jabs_message(MSG_INFO, stderr, "%s", t->help_text);
            if(strcmp(t->name, "help") == 0) {
                size_t i = 0;
                for(const struct help_topic *t2 = topics; t2->name != NULL; t2++) {
                    i++;
                    jabs_message(MSG_INFO, stderr, "%18s", t2->name);
                    if(i % 4 == 0) {
                        jabs_message(MSG_INFO, stderr, "\n");
                    }
                }
                jabs_message(MSG_INFO, stderr, "\nTry also \"help commands\" for a list of possible commands. Try running a command without arguments to get a brief help on usage.");
            }
            jabs_message(MSG_INFO, stderr, "\n");
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
    jabs_message(MSG_INFO, stderr, "JaBS version %s.\n", jabs_version());
    if(git_populated()) {
        jabs_message(MSG_INFO, stderr, "This version of JaBS is compiled from a git repository (branch %s%s).\n", git_branch(), git_dirty() ? ", dirty" : "");
        jabs_message(MSG_INFO, stderr, "Git commit %s dated %s.\n", git_commit_sha1(), git_commit_date());
    }
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_help_commands(script_session *s, int argc, char *const *argv) {
    (void) argv;
    if(argc == 0) {
        jabs_message(MSG_INFO, stderr, "The following commands are available:\n");
        script_print_command_tree(stderr, s->commands);
        return SCRIPT_COMMAND_SUCCESS;
    }
    return SCRIPT_COMMAND_NOT_FOUND;
}
#ifdef JABS_PLUGINS

script_command_status script_identify_plugin(struct script_session *s, int argc, char * const *argv) {
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: identify plugin <path>\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    jabs_plugin *plugin = jabs_plugin_open(argv[0]);
    if(!plugin) {
        jabs_message(MSG_ERROR, stderr, "Plugin %s not found or could not be opened.\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    jabs_message(MSG_INFO, stderr, "Plugin identifies as \"%s\" version \"%s\", type of plugin is %s.\n", plugin->name, plugin->version, jabs_plugin_type_string(plugin->type));
    jabs_plugin_close(plugin);
    return 1;
}

script_command_status script_load_reaction_plugin(script_session *s, int argc, char *const *argv) {
    const int argc_orig = argc;
    struct fit_data *fit = s->fit;
    if(argc < 2) {
        jabs_message(MSG_ERROR, stderr, "Usage: load reaction plugin <file> <target> ...\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    const jibal_isotope *target = jibal_isotope_find(s->jibal->isotopes, argv[1], 0, 0);
    if(!target) {
        jabs_message(MSG_ERROR, stderr, "Could not find isotope matching \"%s\".\n", argv[1]);
        return SCRIPT_COMMAND_FAILURE;
    }
    const char *filename = argv[0];
    jabs_plugin *plugin = jabs_plugin_open(filename);
    if(!plugin) {
        jabs_message(MSG_ERROR, stderr, "Could not load plugin from file \"%s\".\n", filename);
        return EXIT_FAILURE;
    }
    reaction *r = reaction_make(fit->sim->beam_isotope, target, REACTION_PLUGIN, JIBAL_CS_NONE);
    if(!r) {
        jabs_message(MSG_ERROR, stderr, "Could not make a new reaction.\n");
        jabs_plugin_close(plugin);
        return SCRIPT_COMMAND_FAILURE;
    }
    r->plugin = plugin;
    argc += 2;
    argv -= 2;
    jabs_plugin_reaction *pr = jabs_plugin_reaction_init(plugin, s->jibal->isotopes, fit->sim->beam_isotope, target, &argc, &argv);
    if(!pr) {
        jabs_message(MSG_ERROR, stderr, "Plugin failed to initialize.\n");
        reaction_free(r);
        jabs_plugin_close(plugin);
        return SCRIPT_COMMAND_FAILURE;

    }
    r->plugin_r = pr;
    r->product = pr->product;
    r->product_nucleus = pr->product_heavy;
    r->E_min = pr->E_min;
    r->E_max = pr->E_max;
    r->filename = strdup(plugin->filename);
    sim_reactions_add_reaction(fit->sim, r);
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
    jabs_message(MSG_INFO, stderr, "%s\n", cwd);
    free(cwd);
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_cd(struct script_session *s, int argc, char *const *argv) {
    (void) s;
    const int argc_orig = argc;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: cd <path>\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(chdir(argv[0])) {
        jabs_message(MSG_ERROR, stderr, "Could not change directory to \"%s\"\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    argc--;
    argv++;
    return argc_orig - argc;
}
