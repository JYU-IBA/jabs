/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021-2022 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#include <stdlib.h>
#include <string.h>
#include "message.h"
#include "sample.h"
#include "spectrum.h"
#include "simulation.h"
#include "fit.h"
#include "options.h"
#include "jabs.h"
#include "generic.h"
#include "script_command.h"

int script_prepare_sim_or_fit(script_session *s) {
    /* TODO: option to disable RBS or ERD */
    /* TODO: reactions from files? Should "sim_reactions_add" be somewhere else? */
    fit_data *fit = s->fit;
    if(!fit->sm) {
        jabs_message(MSG_ERROR, stderr,"No sample has been defined!\n");
        return -1;
    }
    if(!fit->sim->beam_isotope) {
        jabs_message(MSG_ERROR, stderr,"No ion has been defined!\n");
        return -1;
    }
    if(!fit->sim->det || fit->sim->n_det == 0) {
        jabs_message(MSG_ERROR, stderr,"No detector has been defined!\n");
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
        jabs_message(MSG_ERROR, stderr, "Could not make a sample based on model description. This should never happen.\n");
        return -1;
    }
    if(fit->sim->n_reactions == 0) {
        jabs_message(MSG_WARNING, stderr, "No reactions, adding some automatically. Please be aware there are commands called \"reset reactions\" and \"add reactions\".\n");
        if(fit->sim->rbs) {
            sim_reactions_add_auto(fit->sim, fit->sm, REACTION_RBS, sim_cs(fit->sim, REACTION_RBS)); /* TODO: loop over all detectors and add reactions that are possible (one reaction for all detectors) */
        }
        if(fit->sim->erd) {
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

    jibal_gsto_assign_clear_all(fit->jibal->gsto); /* Is it necessary? No. Here? No. Does it clear old stuff? Yes. */
    if(assign_stopping(fit->jibal->gsto, fit->sim)) {
        jabs_message(MSG_ERROR, stderr,"Could not assign stopping or straggling. Failure. Provide more data, check that JIBAL Z2_max is sufficiently large (currently %i) or disable unwanted reactions (e.g. ERD).\n", s->jibal->config->Z_max);
        return -1;
    }
    jibal_gsto_print_assignments(fit->jibal->gsto);
    jibal_gsto_print_files(fit->jibal->gsto, TRUE);
    jabs_message(MSG_VERBOSE, stderr, "Loading stopping data.\n");
    jibal_gsto_load_all(fit->jibal->gsto);
    simulation_print(stderr, fit->sim);
    s->start = clock();
    return 0;
}

int script_finish_sim_or_fit(script_session *s) {
    s->end = clock();
    double cputime_total = (((double) (s->end - s->start)) / CLOCKS_PER_SEC);
    jabs_message(MSG_INFO, stderr, "...finished! Total CPU time: %.3lf s.\n", cputime_total);

    struct fit_data *fit = s->fit;

    if(fit->sim->n_det == 1) { /* TODO: multidetector automatic spectra saving! */
        size_t i_det = 0;
        sim_workspace *ws = fit_data_ws(fit, i_det);
        if(ws) {
            if(s->output_filename) {
                if(print_spectra(s->output_filename, ws, fit_data_exp(fit, i_det))) {
                    jabs_message(MSG_ERROR, stderr, "Could not save spectra of detector %zu to file \"%s\"\n", i_det, s->output_filename);
                    return EXIT_FAILURE;
                }
            }
        }
    }
    return 0;
}


void script_command_not_found(const char *cmd, const script_command *c) {
    if(c) {
        if(cmd) {
            jabs_message(MSG_ERROR, stderr, "Sub-command \"%s\" is invalid!\n\n", cmd);
        } else {
            jabs_message(MSG_ERROR, stderr, "Not enough arguments!\n\n");
        }
        jabs_message(MSG_ERROR, stderr, "Following subcommands are recognized:\n");
        script_print_commands(stderr, c);
    } else {
        if(cmd) {
            jabs_message(MSG_ERROR, stderr, "Invalid command or argument: \"%s\".\n", cmd);
        } else {
            jabs_message(MSG_ERROR, stderr, "What?\n", cmd);
        }
    }
}

script_command_status script_simulate(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    (void) argc; /* Unused */
    (void) argv; /* Unused */
    if(script_prepare_sim_or_fit(s)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    if(fit_data_workspaces_init(fit)) {
        jabs_message(MSG_ERROR, stderr, "Could not initialize simulation workspace(s).\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    for(size_t i_det = 0; i_det < fit->sim->n_det; i_det++) {
        if(simulate_with_ds(fit->ws[i_det])) {
            jabs_message(MSG_ERROR, stderr, "Simulation failed.\n");
            return SCRIPT_COMMAND_FAILURE;
        }
    }
    script_finish_sim_or_fit(s);
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_fit(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit_data = s->fit;
    if(argc != 1) {
        fprintf(stderr, "Usage: fit [fitvar1,fitvar2,...]\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    fit_params_free(fit_data->fit_params);
    fit_data->fit_params = fit_params_new();
    if(fit_params_add(fit_data->sim, fit_data->sm, fit_data->fit_params, argv[0])) {
        jabs_message(MSG_ERROR, stderr, "Could not add some fit parameters.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(fit_data->fit_params->n == 0) {
        jabs_message(MSG_ERROR, stderr, "No parameters for fit.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
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
    if(fit(fit_data)) {
        jabs_message(MSG_ERROR, stderr, "Fit failed!\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    script_finish_sim_or_fit(s);
    jabs_message(MSG_INFO, stderr, "\nFinal parameters:\n");
    simulation_print(stderr, fit_data->sim);
    jabs_message(MSG_INFO, stderr, "\nFinal profile:\n");
    sample_print(NULL, fit_data->sim->sample, FALSE);
    sample_areal_densities_print(stderr, fit_data->sim->sample, FALSE);
    jabs_message(MSG_INFO, stderr, "\nFinal sample model:\n");
    sample_model_print(NULL, fit_data->sm);
    jabs_message(MSG_INFO, stderr, "\n");
    fit_stats_print(stderr, &fit_data->stats);
    return 0;
}

script_command_status script_save_bricks(script_session *s, int argc, char * const *argv) {
    size_t i_det = 0;
    struct fit_data *fit = s->fit;
    if(script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det) || argc != 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: save bricks [detector] file\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(print_bricks(argv[0], fit_data_ws(fit, i_det))) {
        jabs_message(MSG_ERROR, stderr, "Could not save bricks of detector %zu to file \"%s\"! There should be %zu detector(s).\n", i_det + 1, argv[0], fit->sim->n_det);
        return SCRIPT_COMMAND_FAILURE;
    }
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_save_spectra(script_session *s, int argc, char * const *argv) {
    size_t i_det = 0;
    struct fit_data *fit = s->fit;
    if(script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det) || argc != 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: save spectra [detector] file\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(print_spectra(argv[0], fit_data_ws(fit, i_det), fit_data_exp(fit, i_det))) {
        jabs_message(MSG_ERROR, stderr, "Could not save spectra of detector %zu to file \"%s\"! There should be %zu detector(s).\n", i_det + 1, argv[0], fit->sim->n_det);
        return SCRIPT_COMMAND_FAILURE;
    }
    return SCRIPT_COMMAND_SUCCESS;
}



script_command_status script_save_sample(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit_data = s->fit;
    if(argc != 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: save sample [file]\n");
    }
    if(!fit_data->sm) {
        jabs_message(MSG_ERROR, stderr, "No sample set.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(sample_model_print(argv[0], fit_data->sm)) {
        jabs_message(MSG_ERROR, stderr, "Could not write sample to file \"%s\".\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_save_detector(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit_data = s->fit;
    size_t i_det = 0;
    if(script_get_detector_number(fit_data->sim, TRUE, &argc, &argv, &i_det) || argc != 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: save detector [detector] file\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(detector_print(argv[0], sim_det(fit_data->sim, i_det))) {
        jabs_message(MSG_ERROR, stderr,"Could not write detector %zu to file \"%s\".\n", i_det, argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_remove_reaction(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit_data = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: remove reaction [TYPE] [target_isotope]   OR   remove reaction [number]\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    char *end;
    size_t i = strtoull(argv[0], &end, 10);
    if(*end == '\0') {
        return sim_reactions_remove_reaction(fit_data->sim, i - 1);
    }
    if(argc != 2) {
        jabs_message(MSG_ERROR, stderr, "Usage: remove reaction [TYPE] [target_isotope]   OR   remove reaction [number]\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    reaction_type type = reaction_type_from_string(argv[0]);
    const jibal_isotope *target = jibal_isotope_find(fit_data->jibal->isotopes, argv[1], 0, 0);
    if(type == REACTION_NONE) {
        jabs_message(MSG_ERROR, stderr, "This is not a valid reaction type: \"%s\".\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    if(!target) {
        jabs_message(MSG_ERROR, stderr, "This is not a valid isotope: \"%s\".\n", argv[1]);
        return SCRIPT_COMMAND_FAILURE;
    }
    for(size_t i = 0; i < fit_data->sim->n_reactions; i++) {
        if(fit_data->sim->reactions[i]->type == type && fit_data->sim->reactions[i]->target == target) {
            return sim_reactions_remove_reaction(fit_data->sim, i);
        }
    }
    jabs_message(MSG_ERROR, stderr, "No matching reaction found!\n");
    return SCRIPT_COMMAND_FAILURE;
}

script_command_status script_roi(script_session *s, int argc, char * const *argv) {
    struct roi r;
    if(argc == 3) {
        size_t i = strtoul(argv[0], NULL, 10);
        if(i == 0) {
            jabs_message(MSG_ERROR, stderr, "Detector number must be > 0\n");
            return SCRIPT_COMMAND_FAILURE;
        }
        if(i > s->fit->sim->n_det) {
            jabs_message(MSG_WARNING, stderr, "Warning: Detector number %zu > %zu.\n", i, s->fit->sim->n_det);
        }
        r.i_det = i - 1;
        argv++;
        argc--;
    }
    if (argc == 2) {
        r.i_det = 0;
        r.low = strtoul(argv[0], NULL, 10);
        r.high = strtoul(argv[1], NULL, 10);
    } else {
        jabs_message(MSG_ERROR, stderr, "Usage: roi [det] low high\n");
    }
    fit_data_roi_print(stderr, s->fit, &r);
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_disable(script_session *s, int argc, char *const *argv) {
    if(argc != 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: disable (variable)\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    return script_set_boolean(s, argv[0], FALSE);
}

script_command_status script_enable(script_session *s, int argc, char * const *argv) {
    if(argc != 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: enable (variable)\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    return script_set_boolean(s, argv[0], TRUE);
}

script_command_status script_exit(script_session *s, int argc, char * const *argv) {
    (void) s;
    (void) argc;
    (void) argv;
    return SCRIPT_COMMAND_EXIT;
}

const script_command *script_command_find(const script_command *commands, const char *cmd_string) {
    int found = 0;
    const script_command *c_found = NULL;
    for(const struct script_command *c = commands; c->name != NULL; c++) {
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
    if(found > 1) {
        jabs_message(MSG_ERROR, stderr, "\"%s\" is ambiguous (%i matches):", cmd_string, found);
        for(const struct script_command *c = commands; c->name != NULL; c++) {
            if(strncmp(c->name, cmd_string, strlen(cmd_string)) == 0) {
                jabs_message(MSG_ERROR, stderr, " %s", c->name);
            }
        }
        jabs_message(MSG_ERROR, stderr, "\n");
        return NULL;
    }
    return NULL;
}

script_command_status script_set_boolean(script_session *s, const char *variable, int value) {
    for(jibal_config_var *var = s->cf->vars; var->type != 0; var++) {
        if(strcmp(var->name, variable) == 0) {
            if(var->type != JIBAL_CONFIG_VAR_BOOL) {
                jabs_message(MSG_ERROR, stderr, "Variable exists, but is not boolean.\n");
                return SCRIPT_COMMAND_FAILURE;
            }
            *((int *)var->variable) = value;
            return SCRIPT_COMMAND_SUCCESS;
        }
    }
    return SCRIPT_COMMAND_NOT_FOUND;
}

script_command_status script_execute_command(script_session *s, const char *cmd) {
    int argc = 0;
    script_command_status status;
    char **argv = string_to_argv(cmd, &argc);
    if(!argv) {
        jabs_message(MSG_ERROR, stderr, "Something went wrong in parsing arguments.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc) {
        status = script_execute_command_argv(s, script_commands, argc, argv); /* Note that s->file_depth may be altered (e.g. by script_load_script() */
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
    while(1) {
        if(argc && cmds) { /* Arguments and subcommands remain. Try to find the right one, if possible. */
            const script_command *c = script_command_find(cmds, argv[0]);
            if(c) { /* Subcommand found */
#ifdef DEBUG
                fprintf(stderr, "Debug: Found command %s.\n", c->name);
#endif
                if(c->f) {
#ifdef DEBUG
                    fprintf(stderr, "There is a function %p in command %s. Calling it with %i arguments.\n", (void *) c->f, c->name, argc);
#endif
                    script_command_status status = c->f(s, argc - 1, argv + 1);
                    if(status != SCRIPT_COMMAND_NOT_FOUND) {
                        return status;
                    }
                }
                if(!c->subcommands) {
#ifdef DEBUG
                    fprintf(stderr, "Debug: there are no subcommands.\n");
#endif
                    script_command_not_found(argv[1], NULL);
                    return SCRIPT_COMMAND_NOT_FOUND;
                } else {
                    cmds = c->subcommands;
                    argc--;
                    argv++;
                    c_parent = c;
                    continue;
                }
            } else {
#ifdef DEBUG
                fprintf(stderr, "Debug: Didn't find command %s.\n", argv[0]);
#endif
                script_command_not_found(argv[0], c_parent?c_parent->subcommands:NULL);
            }

        } else {
#ifdef DEBUG
            fprintf(stderr, "Debug: Didn't find command %s (and no more arguments remain).\n", argv[0]);
#endif
            script_command_not_found(NULL, c_parent?c_parent->subcommands:NULL);
        }
        return SCRIPT_COMMAND_NOT_FOUND;
    }
}

script_command_option *options_from_commands(const script_command *commands) {
    size_t n = script_commands_size(commands);
    script_command_option *out = malloc((n+1)*sizeof(script_command_option));
    for(size_t i = 0; i < n; i++) {
        out[i].name = commands[i].name;
        out[i].val = 0;
        out[i].c = &commands[i];
        out[i].var = NULL;
    }
    out[n].name = NULL;
    return out;
}

script_command_option *options_from_jibal_config(jibal_config_var *vars) {
    if(!vars)
        return NULL;
    size_t n = 0;
    for(const jibal_config_var *var = vars; var->type != 0; var++) {
        n++;
    }
    script_command_option *out = malloc((n+1)*sizeof(script_command_option));
    for(size_t i = 0; i < n; i++) {
        out[i].name = vars[i].name;
        out[i].var = &(vars[i]);
        out[i].val = 0;
        out[i].c = NULL;
    }
    out[n].name = NULL;
    return out;
}

script_command_option *options_append(script_command_option *opt_to, const script_command_option *opt_from) { /* Keeps opt_to mostly intact and appends a shallow copy of opt_from to it */
    if(!opt_from)
        return opt_to;
    size_t n_to = options_size(opt_to);
    size_t n_from = options_size(opt_from);
    size_t n = n_to + n_from; /* n real entries */
    opt_to = realloc(opt_to, (n+1) * sizeof(script_command_option)); /* n + 1 allocated (null terminated) */
    for(size_t i = n_to; i < n; i++) {
        opt_to[i] = opt_from[i - n_to];
    }
    opt_to[n].name = NULL;
    return opt_to;
}

void options_sort(script_command_option *opt) {
    qsort((void *)opt, options_size(opt), sizeof(script_command_option), &option_compare);

}


int option_compare(const void *a, const void *b) {
    const script_command_option *o_a = (const script_command_option *)a;
    const script_command_option *o_b = (const script_command_option *)b;
    return strcmp(o_a->name, o_b->name);
}

size_t options_size(const script_command_option *opt) {
    size_t n = 0;
    if(!opt)
        return 0;
    for(const script_command_option *o = opt; o->name != NULL; o++) {
        n++;
    }
    return n;
}

void options_print(FILE *f, const script_command_option *opt) {
    for(const script_command_option *o = opt; o->name != NULL; o++) {
        jabs_message(MSG_INFO, f, "%24s %6i %16p %16p\n", o->name, o->val, o->c, o->var);
    }
}

int script_getopt(script_session *s, int argc, char * const *argv, const script_command_option *options) {
    int found = 0;
    if(!options)
        return -1;
    if(argc == 0)
        return -1;
    const char *a = (argv)[0];
    unsigned long len = strlen(a);
    const script_command_option *opt_found = NULL;
    for(const script_command_option *o = options; o->name != NULL; o++) {
        if(strncmp(o->name, a, len) == 0) { /* At least partial match */
            found++;
            opt_found = o;
            if(strcmp(o->name, a) == 0) { /* Exact match */
                break;
            }
        }
    }
    if(found > 2) {
        jabs_message(MSG_ERROR, stderr, "Ambiguous: %s\n", a);
        return -2;
    }
    if(found == 0) {
        jabs_message(MSG_ERROR, stderr, "No match: %s\n", a);
        return -2;
    }
    if(!opt_found) {
        jabs_message(MSG_ERROR, stderr, "Impossible has happened in script_getopt()\n");
        return -2;
    }
    argc--;
    argv++;
    if(opt_found->c) {
        if(opt_found->c->f) {
            opt_found->c->f(s, argc, argv);
        } else {
            jabs_message(MSG_WARNING, stderr, "Command found, but no there is no function in it.\n");
        }
    } else if(opt_found->var) {
        if(argc < 1) {
            jabs_message(MSG_WARNING, stderr, "Not enough arguments for setting a variable!\n");
        }
        jibal_config_var_set(s->cf->units, opt_found->var, (argv)[0], s->cf->filename);
        argc--;
        argv++;
    }
    return opt_found->val;
}

void script_print_commands(FILE *f, const struct script_command *commands) {
    if(!commands)
        return;
    for(const struct script_command *c = commands; c->name != NULL; c++) {
        if(!c->help_text)
            continue;
        jabs_message(MSG_INFO, f, " %16s    %s\n", c->name, c->help_text);
    }
}

size_t script_commands_size(const script_command *commands) {
    if(!commands)
        return 0;
    size_t n = 0;
    for(const struct script_command *c = commands; c->name != NULL; c++) {
        n++;
    }
    return n;
}

void script_print_command_tree(FILE *f, const struct script_command *commands) {
    const struct script_command *stack[COMMAND_DEPTH];
    const struct script_command *c;
    stack[0] = commands;
    size_t i = 0;
    c = stack[0];
    while(c->name != NULL) {
        while(c->name != NULL) {
            if(c->f) {
                for(size_t j = 1; j <= i; j++) {
                    jabs_message(MSG_INFO, f, "%s ", stack[j]->name);
                }
                jabs_message(MSG_INFO, f, "%s\n", c->name);
            }
            if(c->subcommands && i < (COMMAND_DEPTH - 1)) {
                i++;
                stack[i] = c;
                c = c->subcommands;
            } else {
                c++;
            }
        }
        if(i == 0)
            break;
        c = stack[i];
        i--;
        c++;
    }
}

script_command_status script_load_script(script_session *s, int argc, char * const *argv) {
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: load script [file]\n");
        return EXIT_FAILURE;
    }
    return script_session_load_script(s, argv[0]);
}

script_command_status script_load_sample(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: load sample [file]\n");
        return EXIT_FAILURE;
    }
    sample_model *sm = sample_model_from_file(fit->jibal, argv[1]);
    if(!sm) {
        jabs_message(MSG_INFO, stderr, "Sample load from \"%s\" failed.\n", argv[1]);
        return -1;
    }
    sample_model_free(fit->sm);
    fit->sm = sm;
    return 0;
}

script_command_status script_load_detector(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    size_t i_det = 0;
    if(script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc < 1) {
        jabs_message(MSG_ERROR,  stderr,"Usage: load detector [detector] filename\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    return sim_det_set(fit->sim, detector_from_file(fit->jibal, argv[0]), i_det);
}

script_command_status script_load_experimental(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    size_t i_det = 0;
    if(script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc < 1) {
        jabs_message(MSG_ERROR,  stderr,"Usage: load experimental [detector] filename\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    gsl_histogram *h = spectrum_read(argv[0], sim_det(fit->sim, i_det));
    if(!h) {
        jabs_message(MSG_ERROR,  stderr,"Reading spectrum from file \"%s\" was not successful.\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    if(fit->exp[i_det]) {
        gsl_histogram_free(fit->exp[i_det]);
    }
    fit->exp[i_det] = h;
    return EXIT_SUCCESS;
}

script_command_status script_load_reaction(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: load reaction [file]\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    return sim_reactions_add_r33(fit->sim, fit->jibal->isotopes, argv[0]);
}

script_command_status script_reset_reactions(script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    sim_reactions_free(s->fit->sim);
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_reset_detectors(script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    for(size_t i_det = 0; i_det < s->fit->sim->n_det; i_det++) {
        detector_free(sim_det(s->fit->sim, i_det));
    }
    s->fit->sim->n_det = 0;
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_reset_fit_ranges(script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    fit_data_fit_ranges_free(s->fit);
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_reset_sample(script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    struct fit_data *fit = s->fit;
    sample_model_free(fit->sm);
    fit->sm = NULL;
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_reset_experimental(script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    struct fit_data *fit = s->fit;
    fit_data_exp_free(s->fit);
    fit->exp = calloc(fit->sim->n_det, sizeof(gsl_histogram *));
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_reset(script_session *s, int argc, char * const *argv) {
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
    fit->sim = sim_init(s->jibal);
    fit->exp = calloc(fit->sim->n_det, sizeof(gsl_histogram *));
    jibal_gsto_assign_clear_all(fit->jibal->gsto);
    script_session_reset_vars(s);
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_show_sample(script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    struct fit_data *fit = s->fit;
    if(!fit->sm) {
        jabs_message(MSG_WARNING, stderr, "No sample has been set.\n");
        return SCRIPT_COMMAND_SUCCESS;
    } else {
        return sample_model_print(NULL, fit->sm);
    }
}

script_command_status script_show_simulation(script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    simulation_print(stderr, s->fit->sim);
    return SCRIPT_COMMAND_SUCCESS;
}
script_command_status script_show_fit(script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    fit_data_print(stderr, s->fit);
    return SCRIPT_COMMAND_SUCCESS;
}
script_command_status script_show_detector(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    size_t i_det = 0;
    if(script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    detector_print(NULL, fit->sim->det[i_det]);
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_show_reactions(script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    struct fit_data *fit = s->fit;
    if(fit->sim->n_reactions == 0) {
        jabs_message(MSG_INFO, stderr, "No reactions.\n");
    }
    reactions_print(stderr, fit->sim->reactions, fit->sim->n_reactions);
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_show_variables(script_session *s, int argc, char *const *argv) {
    (void) argc;
    (void) argv;
    const jibal_config_var *var;
    for(var = s->cf->vars; var->type != 0; var++) {
        if(var->variable == NULL)
            continue;
        switch(var->type) {
            case JIBAL_CONFIG_VAR_NONE:
                break;
            case JIBAL_CONFIG_VAR_PATH:
            case JIBAL_CONFIG_VAR_STRING:
                if(*((void **) var->variable) == NULL)
                    continue;
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
    }
    return SCRIPT_COMMAND_SUCCESS;
}
script_command_status script_set(script_session *s, int argc, char * const *argv) {
    script_command_option *opt = NULL;
    opt = options_append(opt, options_from_commands(script_set_commands));
    opt = options_append(opt, options_from_jibal_config(s->cf->vars));
    options_sort(opt); /* TODO: this does not need to be built every time this function is called */
    options_print(stderr, opt);

    if(argc == 0) {
        script_command_not_found(NULL, script_set_commands);
        jabs_message(MSG_INFO, stderr, "\nAlso see 'help set' and 'show variables' for additional variables you can set.\nExample: 'set alpha \"10 deg\"'\n");
    }
    script_getopt(s, argc, argv, opt);
    free(opt);
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_set_ion(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc != 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: set ion [ion]\nExample: set ion 4He\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    const jibal_isotope *isotope = jibal_isotope_find(fit->jibal->isotopes, argv[0], 0, 0);
    if(!isotope) {
        jabs_message(MSG_ERROR, stderr,"No such isotope: %s\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    s->fit->sim->beam_isotope = isotope;
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_set_aperture(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    aperture *a = &fit->sim->beam_aperture;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: set aperture (type) [width|height|diameter (value)] ...\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    while(argc >= 1) {
        for(const jibal_option *o = aperture_option; o->s; o++) { /* look for a matching aperture_option keyword (circle, rectangle) */
            if(strcmp(o->s, argv[0]) == 0) {
                a->type = o->val;
                argv++;
                argc--;
                break;
            }
        }
        if(argc < 2) /* Following things need two arguments */
            break;
        if(strcmp(argv[0], "width") == 0) {
            a->width = jibal_get_val(fit->jibal->units, UNIT_TYPE_DISTANCE, argv[1]);
        } else if(strcmp(argv[0], "height") == 0) {
            a->height = jibal_get_val(fit->jibal->units, UNIT_TYPE_DISTANCE, argv[1]);
        } else if(strcmp(argv[0], "diameter") == 0) {
            a->diameter = jibal_get_val(fit->jibal->units, UNIT_TYPE_DISTANCE, argv[1]);
        } else {
            jabs_message(MSG_ERROR, stderr, "Unrecognized argument (%s)\n", argv[0]);
            return SCRIPT_COMMAND_FAILURE;
        }
        argc -= 2;
        argv += 2;
    }
    if(argc) {
        jabs_message(MSG_ERROR, stderr, "Unexpected extra arguments (%i), starting with %s\n", argc, argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_set_beam(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: set beam energy|fluence|ion (value)...\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    while(argc >= 2) {
        if(strcmp(argv[0], "energy") == 0) {
            fit->sim->beam_E = jibal_get_val(fit->jibal->units, UNIT_TYPE_ENERGY, argv[1]);
        } else if(strcmp(argv[0], "energy_broad") == 0) {
            fit->sim->beam_E_broad = jibal_get_val(fit->jibal->units, UNIT_TYPE_ENERGY, argv[1]);
        } else if(strcmp(argv[0], "fluence") == 0) {
            fit->sim->fluence = jibal_get_val(fit->jibal->units, UNIT_TYPE_ANY, argv[1]);
        } else if(strcmp(argv[0], "ion") == 0) { /* Alternative to 'set ion' is 'set beam ion' here. */
            const jibal_isotope *isotope = jibal_isotope_find(fit->jibal->isotopes, argv[1], 0, 0);
            if(!isotope) {
                jabs_message(MSG_ERROR, stderr, "\"%s\" is not a valid isotope!\n");
                return SCRIPT_COMMAND_FAILURE;
            }
            fit->sim->beam_isotope = isotope;
        }
        argc -= 2;
        argv += 2;
    }
    if(argc) {
        jabs_message(MSG_ERROR, stderr, "Unexpected extra arguments (%i), starting with %s\n", argc, argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_set_detector(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    size_t i_det = 0;
    if(script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc < 2) {
        jabs_message(MSG_ERROR, stderr, "Usage: set detector [number] variable value variable2 value2 ...\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    while(argc >= 2) { /* TODO: replace by argument vector parsing */
        if(detector_set_var(s->jibal, sim_det(fit->sim, i_det), argv[0], argv[1])) {
            jabs_message(MSG_ERROR, stderr, "Can't set \"%s\" to be \"%s\"!\n", argv[0], argv[1]);
            return SCRIPT_COMMAND_FAILURE;
        }
        argv += 2;
        argc -= 2;
    }
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_set_sample(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 2) {
        jabs_message(MSG_ERROR, stderr, "Usage: set sample [sample]\nExample: set sample TiO2 1000tfu Si 10000tfu\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    sample_model *sm_new = sample_model_from_argv(fit->jibal, &argc, &argv);
    if(argc != 0) {
        jabs_message(MSG_ERROR, stderr, "Not all arguments were parsed correctly (starting from \"%s\")! Sample not set.\n", argv[0]);
        sample_model_free(sm_new);
        return SCRIPT_COMMAND_FAILURE;
    }
    sample_model_free(fit->sm);
    fit->sm = sm_new;
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_set_variable(script_session *s, int argc, char * const *argv) {
    if(argc < 1)
        return SCRIPT_COMMAND_FAILURE;
    for(jibal_config_var *var = s->cf->vars; var->type != 0; var++) {
        if(strcmp(var->name, argv[0]) == 0) {
#ifdef DEBUG
            fprintf(stderr, "Setting variable %s to \"%s\"\n", var->name, argv[1]);
#endif
            if(argc != 2) { /* Yes, this check is quite late, but it allows us to check if the variable exists */
                jabs_message(MSG_ERROR, stderr, "Usage: set variable %s [value]\n", var->name);
                return SCRIPT_COMMAND_FAILURE;
            }
            return jibal_config_var_set(s->cf->units, var, argv[1], s->cf->filename);
        }
    }
    return SCRIPT_COMMAND_NOT_FOUND;
}


script_command_status script_add_reaction(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 2) {
        jabs_message(MSG_ERROR, stderr, "Usage: add reaction TYPE isotope\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    reaction *r = reaction_make_from_argv(fit->jibal, fit->sim->beam_isotope, &argc, &argv);
    if(!r) {
        jabs_message(MSG_ERROR, stderr, "Could not make a reaction based on given description.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc != 0) {
        jabs_message(MSG_ERROR, stderr, "Unexpected extra arguments (%i), starting with %s\n", argc, argv[0]);
        reaction_free(r);
        return SCRIPT_COMMAND_FAILURE;
    }
    if(r->cs == JIBAL_CS_NONE) {
        jibal_cross_section_type cs = sim_cs(fit->sim, r->type);
        r->cs = cs;
        jabs_message(MSG_VERBOSE, stderr, "Reaction cross section not given or not valid, assuming default for %s: %s.\n",
                     reaction_name(r), jibal_cs_types[cs].s);
    }
    return sim_reactions_add_reaction(fit->sim, r);
}
script_command_status script_add_reactions(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(!fit->sm) {
        jabs_message(MSG_ERROR, stderr, "Cannot add reactions before sample has been set (I need to know which reactions to add!).\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc > 0) {
        if(strcmp(argv[0], "RBS") == 0) {
            sim_reactions_add_auto(fit->sim, fit->sm, REACTION_RBS, sim_cs(fit->sim, REACTION_RBS));
        } else if(strcmp(argv[0], "ERD") == 0) {
            sim_reactions_add_auto(fit->sim, fit->sm, REACTION_ERD, sim_cs(fit->sim, REACTION_ERD));
        } else {
            jabs_message(MSG_ERROR, stderr, "What is %s anyways? Did you mean \"RBS\" or \"ERD\"?\n", argv[0]);
            return SCRIPT_COMMAND_FAILURE;
        }
        return SCRIPT_COMMAND_SUCCESS;
    }
    if(fit->sim->rbs) {
        sim_reactions_add_auto(fit->sim, fit->sm, REACTION_RBS, sim_cs(fit->sim, REACTION_RBS)); /* TODO: loop over all detectors and add reactions that are possible (one reaction for all detectors) */
    }
    if(fit->sim->erd) {
        sim_reactions_add_auto(fit->sim, fit->sm, REACTION_ERD, sim_cs(fit->sim, REACTION_ERD));
    }
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_add_detector(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR,  stderr,"Usage: add detector filename\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    return fit_data_add_det(fit, detector_from_file(fit->jibal, argv[0])); /* Adds a new detector (and space for experimental spectrum) */
}

script_command_status script_add_fit_range(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    roi range = {.i_det = 0};
    if(argc == 3) { /* Three arguments, first is the detector number */
        range.i_det = strtoul(argv[0], NULL, 10);
        if(range.i_det == 0 || range.i_det > fit->sim->n_det) {
            jabs_message(MSG_ERROR, stderr, "Detector number %zu is not valid (n_det = %zu).\n", range.i_det, fit->sim->n_det);
            return SCRIPT_COMMAND_FAILURE;
        }
        range.i_det--;
        argc--;
        argv++;
    }
    if(argc == 2) {
        range.low = strtoul(argv[0], NULL, 10);
        range.high = strtoul(argv[1], NULL, 10);
    } else {
        jabs_message(MSG_ERROR, stderr, "Usage: add fit_range [detector] low high\n");
        return -1;
    }
    fit_data_fit_range_add(fit, &range);
    return 0;
}

script_command_status script_help(script_session *s, int argc, char * const *argv) {
    (void) fit; /* Unused */
    static const struct help_topic topics[] = {
            {"help", "This is help on help. How meta.\nHelp is available on following topics:\n"},
            {"commands", "I recognize the following commands (try 'help' followed by command name):\n"},
            {"command_tree", "Almost full list of recognized commands:\n"},
            {"version", "JaBS version: "},
            {"set", "The following variables can be set (unit optional, SI units assumed otherwise):\n"},
            {"show", "Show things. Possible things: sim, fit, sample, detector, variables.\n"},
            {"fit", "Make a fit. Provide list of variables to fit.\n"},
            {NULL, NULL}
    };
    if(argc == 0) {
        jabs_message(MSG_INFO, stderr, "Type help [topic] for information on a particular topic or \"help help\" for help on help.\n");
        return 0;
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
                jabs_message(MSG_INFO, stderr, "\n");
            } else if(strcmp(t->name, "commands") == 0) {
                script_print_commands(stderr, script_commands);
            } else if(strcmp(t->name, "command_tree") == 0) {
                script_print_command_tree(stderr, script_commands);
            } else if(strcmp(t->name, "version") == 0) {
                jabs_message(MSG_INFO, stderr, "%s\n", jabs_version());
            } else if(strcmp(t->name, "set") == 0) {
                if(!s->cf || !s->cf->vars)
                    break;
                size_t i = 0;
                for(jibal_config_var *var = s->cf->vars; var->type != JIBAL_CONFIG_VAR_NONE; var++) {
                    if(var->type != JIBAL_CONFIG_VAR_UNIT)
                        continue;
                    i++;
                    jabs_message(MSG_INFO, stderr," %25s", var->name);
                    if(i % 3 == 0) {
                        jabs_message(MSG_INFO, stderr,"\n");
                    }
                }
                fprintf(stderr, "\n\nThe following variables are not in SI units:\n");
                for(jibal_config_var *var = s->cf->vars; var->type != JIBAL_CONFIG_VAR_NONE; var++) {
                    if(var->type == JIBAL_CONFIG_VAR_UNIT)
                        continue;
                    jabs_message(MSG_INFO, stderr, " %25s: %s", var->name, jibal_config_var_type_name(var->type));
                    if(var->type == JIBAL_CONFIG_VAR_OPTION && var->option_list) {
                        jabs_message(MSG_INFO, stderr, " (");
                        for(const jibal_option *o = var->option_list; o->s; o++) {
                            jabs_message(MSG_INFO, stderr, "%s%s", o == var->option_list ? "":", ", o->s);
                        }
                        jabs_message(MSG_INFO, stderr, ")\n");
                    } else {
                        jabs_message(MSG_INFO, stderr, "\n");
                    }
                }
                jabs_message(MSG_INFO, stderr,"\n\nAlso the following things can be set (special syntax applies for each):\n");
                script_print_commands(stderr, script_set_commands);
            }
            return 0;
        }
    }

    for(const struct script_command *c = script_commands; c->name != NULL; c++) {
        if(strcmp(c->name, argv[0]) == 0) {
            if(!found) { /* There wasn't a help topic  */
                jabs_message(MSG_INFO, stderr, "\"%s\" is a valid command, but no additional help is available!\n\n", c->name);
            }
            found++;
            if(c->subcommands) {
                jabs_message(MSG_INFO, stderr, "At least the following sub-commands are recognized:\n", c->name);
                script_print_commands(stderr, c->subcommands);
            }
            break;
        }
    }

    if(!found) {
        jabs_message(MSG_ERROR, stderr,"Sorry, no help for '%s' available.\n", argv[0]);
        return SCRIPT_COMMAND_NOT_FOUND;
    }
    return SCRIPT_COMMAND_SUCCESS;
}


