/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021-2022 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_SCRIPT_COMMAND_H
#define JABS_SCRIPT_COMMAND_H


struct script_session;

#include "script_session.h"

#define COMMAND_DEPTH 8

typedef enum script_command_status {
    SCRIPT_COMMAND_SUCCESS = EXIT_SUCCESS,
    SCRIPT_COMMAND_FAILURE = EXIT_FAILURE,
    SCRIPT_COMMAND_NOT_FOUND = 101,
    SCRIPT_COMMAND_EXIT = 102,
    SCRIPT_COMMAND_EOF = 103
} script_command_status;

typedef struct script_command {
    const char *name;
    script_command_status (*f)(struct script_session *, int, char * const *);
    const char *help_text; /* Short help text */
    const struct script_command *subcommands;
    jibal_config_var *var;
    int val;
} script_command;

struct help_topic {
    const char *name;
    const char *help_text;
};

int script_getopt(struct script_session *s, const script_command *commands, int *argc, char *const **argv, script_command_status *status_out); /* Parses argument vector, finds script_command_option "c" and calls "c->f()" or sets c->var. */

struct script_command *script_commands_create(struct script_session *s);
script_command *script_commands_append(script_command *c_to, const script_command *c_from);
void script_commands_sort(script_command *commands);
int command_compare(const void *a, const void *b);
script_command *script_commands_from_jibal_config(jibal_config_var *vars);
void script_commands_print(FILE *f, const struct script_command *commands);
size_t script_commands_size(const script_command *commands);
void script_print_command_tree(FILE *f, const struct script_command *commands);
script_command_status script_execute_command(struct script_session *s, const char *cmd);
script_command_status script_execute_command_argv(struct script_session *s, const script_command *commands, int argc, char **argv);
void script_command_not_found(const char *cmd, const script_command *c);
const script_command *script_command_find(const script_command *commands, const char *cmd_string);
script_command_status script_set_boolean(struct script_session *s, const char *variable, int value);

script_command_status script_add_detector(struct script_session *s, int argc, char * const *argv);
script_command_status script_add_fit_range(struct script_session *s, int argc, char * const *argv);
script_command_status script_add_reaction(struct script_session *s, int argc, char * const *argv);
script_command_status script_add_reactions(struct script_session *s, int argc, char * const *argv);
script_command_status script_disable(struct script_session *s, int argc, char * const *argv);
script_command_status script_enable(struct script_session *s, int argc, char * const *argv);
script_command_status script_exit(struct script_session *s, int argc, char * const *argv);
script_command_status script_fit(struct script_session *s, int argc, char * const *argv);
script_command_status script_help(struct script_session *s, int argc, char * const *argv);
script_command_status script_load_detector(struct script_session *s, int argc, char *const *argv);
script_command_status script_load_experimental(struct script_session *s, int argc, char *const *argv);
script_command_status script_load_reaction(struct script_session *s, int argc, char *const *argv);
script_command_status script_load_sample(struct script_session *s, int argc, char *const *argv);
script_command_status script_load_script(struct script_session *s, int argc, char *const *argv);
script_command_status script_remove_reaction(struct script_session *s, int argc, char * const *argv);
script_command_status script_reset(struct script_session *s, int argc, char * const *argv);
script_command_status script_reset_detectors(struct script_session *s, int argc, char * const *argv);
script_command_status script_reset_experimental(struct script_session *s, int argc, char * const *argv);
script_command_status script_reset_fit_ranges(struct script_session *s, int argc, char * const *argv);
script_command_status script_reset_reactions(struct script_session *s, int argc, char * const *argv);
script_command_status script_reset_sample(struct script_session *s, int argc, char * const *argv);
script_command_status script_roi(struct script_session *s, int argc, char * const *argv);
script_command_status script_save_bricks(struct script_session *s, int argc, char * const *argv);
script_command_status script_save_detector(struct script_session *s, int argc, char * const *argv);
script_command_status script_save_sample(struct script_session *s, int argc, char * const *argv);
script_command_status script_save_spectra(struct script_session *s, int argc, char * const *argv);
script_command_status script_show_detector(struct script_session *s, int argc, char * const *argv);
script_command_status script_show_fit(struct script_session *s, int argc, char * const *argv);
script_command_status script_show_reactions(struct script_session *s, int argc, char * const *argv);
script_command_status script_show_sample(struct script_session *s, int argc, char * const *argv);
script_command_status script_show_simulation(struct script_session *s, int argc, char * const *argv);
script_command_status script_show_variables(struct script_session *s, int argc, char * const *argv);
script_command_status script_simulate(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_aperture(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_beam(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_ion(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_detector(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_sample(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_variable(struct script_session *s, int argc, char * const *argv);

#endif //JABS_SCRIPT_COMMAND_H
