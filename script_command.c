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
        script_commands_print(stderr, c);
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
    if(found > 1) {
        jabs_message(MSG_ERROR, stderr, "\"%s\" is ambiguous (%i matches):", cmd_string, found);
        for(const script_command *c = commands; c; c = c->next) {
            if(strncmp(c->name, cmd_string, strlen(cmd_string)) == 0) {
                jabs_message(MSG_ERROR, stderr, " %s", c->name);
            }
        }
        jabs_message(MSG_ERROR, stderr, "\n");
        return NULL;
    }
    return NULL;
}

script_command_status script_set_var(struct script_session *s, jibal_config_var *var, int argc,  char * const *argv) {
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Not enough arguments to set variable.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    jibal_config_var_set(s->jibal->units, var, argv[0], NULL);
    return 1; /* Number of arguments */
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
    while(1) {
        if(argc && cmds) { /* Arguments and subcommands remain. Try to find the right one, if possible. */
            const script_command *c = script_command_find(cmds, argv[0]);
            if(c) { /* Subcommand found */
#ifdef DEBUG
                fprintf(stderr, "Debug: Found command %s.\n", c->name);
#endif
                if(c->f) {
#ifdef DEBUG
                    fprintf(stderr, "There is a function %p in command %s. Calling it with %i arguments.\n", (void *) c->f, c->name, argc - 1);
#endif
                    script_command_status status = c->f(s, argc - 1, argv + 1);
                    if(status != SCRIPT_COMMAND_NOT_FOUND) {
                        return status;
                    }
                } else if(c->var) {
                    if(!c_parent) {
                        jabs_message(MSG_ERROR, stderr, "Command/option \"%s\" is a variable, but there is no parent command. This is highly unusual.\n", c->name);
                        return SCRIPT_COMMAND_FAILURE;
                    }
                    if(!c_parent->f_var) {
                        jabs_message(MSG_ERROR, stderr, "Command/option \"%s\" is a variable, but there is no function to handle variables in parent command (\"%s\"). This is highly unusual.\n", c->name, c_parent->name);
                        return SCRIPT_COMMAND_FAILURE;
                    }
                    c_parent->f_var(s, c->var, argc - 1, argv + 1);
                }
                if(c->subcommands) {
                    cmds = c->subcommands;
                    argc--;
                    argv++;
                    c_parent = c;
                    continue;
                }  else {
#ifdef DEBUG
                    fprintf(stderr, "Debug: there are no subcommands or a function / variable in \"%s\". There is a val: %i\n", c->name, c->val);
#endif
                    script_command_not_found(argv[1], NULL);
                    return SCRIPT_COMMAND_NOT_FOUND;
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

int command_compare(const void *a, const void *b) {
    const script_command *o_a = (const script_command *)a;
    const script_command *o_b = (const script_command *)b;
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

int script_getopt(script_session *s, const script_command *commands, int *argc, char *const **argv, script_command_status *status_out) {
    int found = 0;
    if(!commands || !argc || !argv) {
        *status_out = SCRIPT_COMMAND_FAILURE;
        return -1;
    }
    if(*argc == 0) { /* Nothing remains */
        *status_out = SCRIPT_COMMAND_SUCCESS;
        return -1;
    }
    const char *a = (*argv)[0];
    size_t len = strlen(a);
    const script_command *c_found = NULL;
    for(const script_command *c = commands; c; c = c->next) {
        if(strncmp(c->name, a, len) == 0) { /* At least partial match */
            found++;
            c_found = c;
            if(strcmp(c->name, a) == 0) { /* Exact match */
                break;
            }
        }
    }
    if(found > 1) {
        jabs_message(MSG_ERROR, stderr, "Ambiguous: %s. Matches:\n", a);
        for(const script_command *c = commands; c; c = c->next) {
            if(strncmp(c->name, a, len) == 0) {
                jabs_message(MSG_ERROR, stderr, " %s\n", c->name);
            }
        }
        jabs_message(MSG_ERROR, stderr, "\n");
        return -1;
    }
    if(found == 0) {
        jabs_message(MSG_ERROR, stderr, "Could not find an option or command with \"%s\".\n", a);
        *status_out = SCRIPT_COMMAND_NOT_FOUND;
        return -1;
    }
#ifdef DEBUG
    fprintf(stderr, "Found exactly 1 hit (%s matched with %s). f = %p, var = %p, val = %i\n", c_found->name, a, (void *)c_found->f, (void *)c_found->var, c_found->val);
#endif
    (*argc)--;
    (*argv)++;
    if(c_found->f) {
        *status_out = c_found->f(s, *argc, *argv);
    } else if(c_found->var) {
        if(*argc < 1) {
            jabs_message(MSG_WARNING, stderr, "Not enough arguments for setting a variable (%s)!\n", c_found->var->name);
            *status_out = SCRIPT_COMMAND_FAILURE;
            return -1;
        }
        jibal_config_var_set(s->cf->units, c_found->var, (*argv)[0], s->cf->filename);
        *status_out = SCRIPT_COMMAND_SUCCESS;
        (*argc)--;
        (*argv)++;
    }
    return c_found->val;
}

script_command *script_command_new(const char *name, const char *help_text, int val, script_command_status (*f)(struct script_session *, int, char * const *)) {
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
int script_command_set_function(script_command *c, script_command_status (*f)(struct script_session *, int, char * const *)) {
    if(!c)
        return EXIT_FAILURE;
    c->f = f;
    if(c->var) {
        free(c->var);
        c->var = NULL;
    }
    return EXIT_SUCCESS;
}

int script_command_set_var(script_command *c, jibal_config_var_type type, const void *variable, const jibal_option *option_list) {
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
#ifdef DEBUG
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
    struct script_command *stack[COMMAND_DEPTH];
    struct script_command *c, *c_old = NULL;
    stack[0] = head;
    size_t i = 0;
    c = stack[0];
    while(c) {
        if(c->subcommands && i < COMMAND_DEPTH) { /* Go deeper, push existing pointer to stack */
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

void script_command_list_add_command(script_command **head, script_command *c_new) {
    if(!c_new)
        return;
    if(*head == NULL) {
        *head = c_new;
        return;
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


    script_command *c = script_command_new("help", "Help.", 0, &script_help);
    script_command_list_add_command(&head, c);

    script_command *c_set = script_command_new("set", "Set something.", 0, NULL);
    c_set->f_var = &script_set_var;
    script_command_list_add_command(&head, c_set);
    script_command_list_add_command(&c_set->subcommands, script_command_new("aperture", "Set aperture.", 0, &script_set_aperture));
    script_command_list_add_command(&c_set->subcommands, script_command_new("beam", "Set beam properties.", 0, &script_set_beam));
#if 0
    script_command_list_add_command(&c_set->subcommands, script_command_new("detector", "Set detector properties.", 0, &script_set_detector));
#endif
    script_command_list_add_command(&c_set->subcommands, script_command_new("ion", "Set incident ion (isotope).", 0, &script_set_ion));
    script_command_list_add_command(&c_set->subcommands, script_command_new("sample", "Set sample.", 0, &script_set_sample));
    script_command_list_add_command(&c_set->subcommands, script_command_new("variable", "Set a variable.", 0, &script_set_variable)); /* TODO: unnecessary / incomplete, remove */

    const jibal_config_var vars[] = {
            {JIBAL_CONFIG_VAR_UNIT,   "fluence",              &sim->fluence,                       NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "energy",               &sim->beam_E,                        NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "energy_broad",         &sim->beam_E_broad,                  NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "emin",                 &sim->emin,                          NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "alpha",                &sim->sample_theta,                  NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "phi",                  &sim->sample_phi,                    NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "channeling",           &sim->channeling_offset,             NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "channeling_slope",     &sim->channeling_slope,              NULL},
            {JIBAL_CONFIG_VAR_STRING, "output",               &s->output_filename,                 NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "erd",                  &sim->erd,                           NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "rbs",                  &sim->rbs,                           NULL},
            {JIBAL_CONFIG_VAR_SIZE,   "fit_maxiter",          &fit->n_iters_max,                   NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "fit_xtol",             &fit->xtol,                          NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "fit_gtol",             &fit->gtol,                          NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "fit_ftol",             &fit->ftol,                          NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "ds",                   &sim->params.ds,                     NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "rk4",                  &sim->params.rk4,                    NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "stop_step_incident",   &sim->params.stop_step_incident,     NULL},
            {JIBAL_CONFIG_VAR_UNIT,   "stop_step_exiting",    &sim->params.stop_step_exiting,      NULL},
            {JIBAL_CONFIG_VAR_DOUBLE, "stop_step_fudge",      &sim->params.stop_step_fudge_factor, NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "nucl_stop_accurate",   &sim->params.nucl_stop_accurate,     NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "mean_conc_and_energy", &sim->params.mean_conc_and_energy,   NULL},
            {JIBAL_CONFIG_VAR_BOOL,   "geostragg",            &sim->params.geostragg,              NULL},
            {JIBAL_CONFIG_VAR_NONE, NULL, NULL,                                                      NULL}
    };
    c = script_command_list_from_vars_array(vars, 0);
    script_command_list_add_command(&c_set->subcommands, c);

    script_command *c_load = script_command_new("load", "Load something.", 0, NULL);
    script_command_list_add_command(&head, c_load);
    script_command_list_add_command(&c_load->subcommands, script_command_new("detector", "Load (replace) a detector.", 0, &script_load_detector));
    script_command_list_add_command(&c_load->subcommands, script_command_new("experimental", "Load an experimental spectrum.", 0, &script_load_experimental));
    script_command_list_add_command(&c_load->subcommands, script_command_new("script", "Load (run) a script.", 0, &script_load_script));
    script_command_list_add_command(&c_load->subcommands, script_command_new("sample", "Load a sample.", 0, &script_load_sample));
    script_command_list_add_command(&c_load->subcommands, script_command_new("reaction",  "Load a reaction from R33 file.", 0, &script_load_reaction));

    script_command *c_show = script_command_new("show", "Show information on things.", 0, NULL);
    script_command_list_add_command(&head, c_show);
    script_command_list_add_command(&c_show->subcommands, script_command_new("detector", "Show detector.", 0, &script_show_detector));
    script_command_list_add_command(&c_show->subcommands, script_command_new("fit" ,"Show fit." , 0, &script_show_fit));
    script_command_list_add_command(&c_show->subcommands, script_command_new("reactions","Show reactions." , 0, &script_show_reactions));
    script_command_list_add_command(&c_show->subcommands, script_command_new("sample", "Show sample.", 0, &script_show_sample));
    script_command_list_add_command(&c_show->subcommands, script_command_new("simulation", "Show simulation.", 0, &script_show_simulation));
    script_command_list_add_command(&c_show->subcommands, script_command_new("variables" , "Show variables.", 0, &script_show_variables));

    script_command_list_add_command(&head, script_command_new("exit", "Exit.", 0, &script_exit));

    c = script_command_new("save", "Save something.", 0, NULL);
    script_command_list_add_command(&head, c);
    script_command_list_add_command(&c->subcommands, script_command_new("bricks", "Save bricks.", 0, &script_save_bricks));
    script_command_list_add_command(&c->subcommands, script_command_new("detector", "Save detector.", 0, &script_save_detector));
    script_command_list_add_command(&c->subcommands, script_command_new("sample", "Save sample.", 0, &script_save_sample));
    script_command_list_add_command(&c->subcommands, script_command_new("spectra", "Save spectra.", 0, &script_save_spectra));

    c = script_command_new("add", "Add something.", 0, NULL);
    script_command_list_add_command(&head, c);
    script_command_list_add_command(&c->subcommands, script_command_new("detector", "Add a detector.", 0, &script_add_detector));
    script_command_list_add_command(&c->subcommands, script_command_new("fit_range", "Add a fit range", 0, &script_add_fit_range));
    script_command_list_add_command(&c->subcommands, script_command_new("reaction", "Add a reaction.", 0, &script_add_reaction));
    script_command_list_add_command(&c->subcommands, script_command_new("reactions", "Add reactions (of some type).", 0, &script_add_reactions));

    c = script_command_new("remove", "Remove something.", 0, NULL);
    script_command_list_add_command(&head, c);
    script_command_list_add_command(&c->subcommands, script_command_new("reaction", "Remove reaction.", 0, &script_remove_reaction));

    c = script_command_new("reset", "Reset something (or everything).", 0, NULL);
    script_command_list_add_command(&head, c);
    script_command_list_add_command(&c->subcommands, script_command_new("detectors", "Reset detectors.", 0, &script_reset_detectors));
    script_command_list_add_command(&c->subcommands, script_command_new("experimental", "Reset experimental spectra.", 0, &script_reset_experimental));
    script_command_list_add_command(&c->subcommands, script_command_new("fit_ranges", "Reset fit ranges.", 0, &script_reset_fit_ranges));
    script_command_list_add_command(&c->subcommands, script_command_new("reactions", "Reset reactions.", 0, &script_reset_reactions));
    script_command_list_add_command(&c->subcommands, script_command_new("sample", "Reset sample.", 0, &script_reset_sample));

    c = script_command_new("fit", "Do a fit.", 0, script_fit);
    script_command_list_add_command(&head, c);

    c = script_command_new("enable", "Set boolean variable to true.", 0, script_enable);
    script_command_list_add_command(&head, c);

    c = script_command_new("disable", "Set boolean variable to false.", 0, script_disable);
    script_command_list_add_command(&head, c);

    c = script_command_new("roi", "Show information from a region of interest.", 0, script_roi);
    script_command_list_add_command(&head, c);

    c = script_command_new("simulate", "Run a simulation.", 0, script_simulate);
    script_command_list_add_command(&head, c);

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
    const struct script_command *stack[COMMAND_DEPTH];
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
        if(c->subcommands && i < COMMAND_DEPTH) { /* Go deeper, push existing pointer to stack */
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

script_command_status script_load_script(script_session *s, int argc, char * const *argv) {
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: load script [file]\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    return script_session_load_script(s, argv[0]);
}

script_command_status script_load_sample(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: load sample [file]\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    sample_model *sm = sample_model_from_file(fit->jibal, argv[1]);
    if(!sm) {
        jabs_message(MSG_INFO, stderr, "Sample load from \"%s\" failed.\n", argv[1]);
        return -1;
    }
    sample_model_free(fit->sm);
    fit->sm = sm;
    if(s->fit->sim->n_reactions > 0) {
        jabs_message(MSG_WARNING, stderr, "Reactions were reset automatically, since the sample was changed.\n");
        sim_reactions_free(fit->sim);
    }
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
    return 1; /* Number of arguments */
}

script_command_status script_load_reaction(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: load reaction [file]\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(sim_reactions_add_r33(fit->sim, fit->jibal->isotopes, argv[0])) {
        return SCRIPT_COMMAND_FAILURE;
    } else {
        return 1;
    }
}

script_command_status script_reset_reactions(script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    sim_reactions_free(s->fit->sim);
    return 0;
}

script_command_status script_reset_detectors(script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    for(size_t i_det = 0; i_det < s->fit->sim->n_det; i_det++) {
        detector_free(sim_det(s->fit->sim, i_det));
    }
    s->fit->sim->n_det = 0;
    return 0;
}

script_command_status script_reset_fit_ranges(script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    fit_data_fit_ranges_free(s->fit);
    return 0;
}

script_command_status script_reset_sample(script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    struct fit_data *fit = s->fit;
    sample_model_free(fit->sm);
    fit->sm = NULL;
    return 0;
}

script_command_status script_reset_experimental(script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    struct fit_data *fit = s->fit;
    fit_data_exp_free(s->fit);
    fit->exp = calloc(fit->sim->n_det, sizeof(gsl_histogram *));
    return 0;
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
    return 0;
}

script_command_status script_show_sample(script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    struct fit_data *fit = s->fit;
    if(!fit->sm) {
        jabs_message(MSG_WARNING, stderr, "No sample has been set.\n");
        return 0;
    } else {
#ifdef DEBUG
        char *sample_str = sample_model_to_string(fit->sm);
        fprintf(stderr, "Sample: %s\n", sample_str);
        free(sample_str);
#endif
        if(sample_model_print(NULL, fit->sm)) {
            return SCRIPT_COMMAND_FAILURE;
        } else {
            return 0;
        }
    }
}

script_command_status script_show_simulation(script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    simulation_print(stderr, s->fit->sim);
    return 0;
}
script_command_status script_show_fit(script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    fit_data_print(stderr, s->fit);
    return 0;
}
script_command_status script_show_detector(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    size_t i_det = 0;
    int argc_orig = argc;
    if(script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    detector_print(NULL, fit->sim->det[i_det]);
    return argc_orig - argc; /* Number of arguments */
}

script_command_status script_show_reactions(script_session *s, int argc, char * const *argv) {
    (void) argc;
    (void) argv;
    struct fit_data *fit = s->fit;
    if(fit->sim->n_reactions == 0) {
        jabs_message(MSG_INFO, stderr, "No reactions.\n");
    }
    reactions_print(stderr, fit->sim->reactions, fit->sim->n_reactions);
    return 0;
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
    return 0;
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
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: set aperture (type) [width|height|diameter (value)] ...\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    aperture *a = aperture_from_argv(s->jibal, &argc, &argv);
    if(!a) {
        jabs_message(MSG_ERROR, stderr, "Aperture could not be parsed.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(argc) {
        free(a);
        jabs_message(MSG_ERROR, stderr, "Unexpected extra arguments (%i), starting with %s\n", argc, argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    fit->sim->beam_aperture = a;
    return SCRIPT_COMMAND_SUCCESS;
}

script_command_status script_set_beam(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    int argc_consumed = 0;
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: set beam energy|fluence|ion (value)...\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    while(argc >= 2) { /* TODO: replace with something */
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
        } else {
            break;
        }
        argc -= 2;
        argv += 2;
        argc_consumed += 2;
    }
    return argc_consumed;
}

#if 0
script_command_status script_set_detector(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    static const script_command extra_commands[] = {
            {"aperture",    NULL, "Aperture description", NULL, NULL, 'a'},
            {"calibration", NULL, "Calibration",          NULL, NULL, 'c'},
            {"foil",        NULL, NULL,                   NULL, NULL, 'f'},
            {NULL, 0,             NULL,                   NULL, NULL, 0}
    };
    size_t i_det = 0;
    if(script_get_detector_number(fit->sim, TRUE, &argc, &argv, &i_det)) {
        return SCRIPT_COMMAND_FAILURE;
    }
    script_command *commands = NULL;
    detector *det = sim_det(fit->sim, i_det);
    jibal_config_var *vars = detector_make_vars(det);
    commands = script_commands_append(commands, script_commands_from_jibal_config(vars));
    commands = script_commands_append(commands, extra_commands);
    script_commands_sort(commands);
    //options_print(stderr, opt);
    if(argc < 1) {
        jabs_message(MSG_ERROR, stderr, "Usage: set detector [number] variable value variable2 value2 ...\n");
        char *matches = script_commands_list_matches(commands, "");
        jabs_message(MSG_ERROR, stderr, "These variables can be set: %s\n", matches);
        free(matches);
        free(vars);
        free(commands);
        return SCRIPT_COMMAND_FAILURE;
    }

    script_command_status status = SCRIPT_COMMAND_SUCCESS;
    while(argc >= 1 && status == SCRIPT_COMMAND_SUCCESS) {
        int val = script_getopt(s, commands, &argc, &argv, &status);
        if(val < 0)
            break;
        switch(val) {
            case 'a':
                if(detector_aperture_set_from_argv(s->jibal, det, &argc, &argv)) {
                    status = SCRIPT_COMMAND_FAILURE;
                }
                break;
            case 'f':
                if(detector_foil_set_from_argv(s->jibal, det, &argc, &argv)) {
                    status = SCRIPT_COMMAND_FAILURE;
                }
                break;
            default:
                break;
        }
    }
    free(commands);
    free(vars);
    if(argc != 0) {
        jabs_message(MSG_ERROR, stderr, "Extra arguments, not enough arguments or generic parse error starting at \"%s\"\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    if(status == SCRIPT_COMMAND_NOT_FOUND)
        status = SCRIPT_COMMAND_FAILURE;
    return status;
}
#endif

script_command_status script_set_sample(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 2) {
        jabs_message(MSG_ERROR, stderr, "Usage: set sample [sample]\nExample: set sample TiO2 1000tfu Si 10000tfu\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    int argc_orig = argc;
    sample_model *sm_new = sample_model_from_argv(fit->jibal, &argc, &argv);
    int argc_consumed = argc_orig - argc;
    sample_model_free(fit->sm);
    fit->sm = sm_new;
    if(s->fit->sim->n_reactions > 0) {
        jabs_message(MSG_WARNING, stderr, "Reactions were reset automatically, since the sample was changed.\n");
        sim_reactions_free(fit->sim);
    }
    return argc_consumed;
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
            if(jibal_config_var_set(s->cf->units, var, argv[1], s->cf->filename)) {
                return SCRIPT_COMMAND_FAILURE;
            } else {
                return 2; /* Number of arguments consumed */
            }
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
    int argc_orig = argc;
    reaction *r = reaction_make_from_argv(fit->jibal, fit->sim->beam_isotope, &argc, &argv);
    int argc_consumed = argc_orig - argc;
    if(!r) {
        jabs_message(MSG_ERROR, stderr, "Could not make a reaction based on given description.\n");
        return SCRIPT_COMMAND_FAILURE;
    }
    if(r->cs == JIBAL_CS_NONE) {
        jibal_cross_section_type cs = sim_cs(fit->sim, r->type);
        r->cs = cs;
        jabs_message(MSG_VERBOSE, stderr, "Reaction cross section not given or not valid, assuming default for %s: %s.\n",
                     reaction_name(r), jibal_cs_types[cs].s);
    }
    if(sim_reactions_add_reaction(fit->sim, r)) {
        return SCRIPT_COMMAND_FAILURE;
    } else {
        return argc_consumed;
    }
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
        return 1;
    }
    if(fit->sim->rbs) {
        sim_reactions_add_auto(fit->sim, fit->sm, REACTION_RBS, sim_cs(fit->sim, REACTION_RBS)); /* TODO: loop over all detectors and add reactions that are possible (one reaction for all detectors) */
    }
    if(fit->sim->erd) {
        sim_reactions_add_auto(fit->sim, fit->sm, REACTION_ERD, sim_cs(fit->sim, REACTION_ERD));
    }
    return 0;
}

script_command_status script_add_detector(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    if(argc < 1) {
        jabs_message(MSG_ERROR,  stderr,"Usage: add detector filename\n", argv[0]);
        return SCRIPT_COMMAND_FAILURE;
    }
    if(fit_data_add_det(fit, detector_from_file(fit->jibal, argv[0]))) { /* Adds a new detector (and space for experimental spectrum) */
        return SCRIPT_COMMAND_FAILURE;
    } else {
        return 1;
    }
}

script_command_status script_add_fit_range(script_session *s, int argc, char * const *argv) {
    struct fit_data *fit = s->fit;
    roi range = {.i_det = 0};
    if(argc == 3) { /* Three arguments, first is the detector number */ /* TODO: we can't rely on fixed number of arguments! */
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
    return argc; /* TODO: we can't rely on fixed number of arguments! */
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
        return SCRIPT_COMMAND_FAILURE;
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
                script_commands_print(stderr, s->commands);
            } else if(strcmp(t->name, "command_tree") == 0) {
                script_print_command_tree(stderr, s->commands);
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
            }
            return 0;
        }
    }

    for(const script_command *c = s->commands; c; c = c->next) { /* TODO: deeper? */
        if(strcmp(c->name, argv[0]) == 0) {
            if(!found) { /* There wasn't a help topic  */
                jabs_message(MSG_INFO, stderr, "\"%s\" is a valid command, but no additional help is available!\n\n", c->name);
            }
            found++;
            if(c->subcommands) {
                jabs_message(MSG_INFO, stderr, "At least the following sub-commands are recognized:\n", c->name);
                script_commands_print(stderr, c->subcommands);
            }
            break;
        }
    }

    if(!found) {
        jabs_message(MSG_ERROR, stderr,"Sorry, no help for '%s' available.\n", argv[0]);
        return SCRIPT_COMMAND_NOT_FOUND;
    }
    return 1; /* TODO: this stuff only works with one argument, right? */
}


