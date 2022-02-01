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

#define SCRIPT_COMMAND_SUCCESS (0) /* Anything above and including zero is success */
#define SCRIPT_COMMAND_FAILURE (-1)
#define SCRIPT_COMMAND_NOT_FOUND (-2)
#define SCRIPT_COMMAND_EXIT (-3)
#define SCRIPT_COMMAND_EOF (-4)

typedef int script_command_status; /* Script commands should return negative on error (see defines above) and number of arguments (zero or positive) consumed successfully */

typedef struct script_command {
    char *name;
    script_command_status (*f)(struct script_session *, int, char * const *); /* Function to process argument vectors. Called if it is non-NULL, even if subcommands exist! Return value is important.*/
    script_command_status (*f_var)(struct script_session *, jibal_config_var *var, int, char * const *); /* Function to process variables (from subcommands). */
    script_command_status (*f_val)(struct script_session *, int, int, char * const *); /* Function to process vals (from subcommands). */
    char *help_text; /* Short help text */
    struct script_command *subcommands;
    jibal_config_var *var;
    int val;
    struct script_command *next;
} script_command;

struct help_topic {
    const char *name;
    const char *help_text;
};


const char *script_command_status_to_string(script_command_status status);

script_command *script_command_new(const char *name, const char *help_text, int val, script_command_status (*f)(struct script_session *, int, char * const *)); /* Allocates new command that doesn't do anything. Can return var. */
int script_command_set_var(script_command *c, jibal_config_var_type type, void *variable, const jibal_option *option_list);
void script_command_free(script_command *c);
void script_commands_free(script_command *head);

script_command *script_command_list_find_tail(script_command *head);
void script_command_list_add_command(script_command **head, script_command *c_new);
script_command *script_command_list_from_command_array(const script_command *commands); /* commands is an array, must be "null-terminated" (name of last element is NULL pointer). Deep copy will be made. */
script_command *script_command_list_from_vars_array(const jibal_config_var *vars, jibal_config_var_type type); /* vars is an array, must be "null-terminated" (name of last element is NULL pointer). Deep copy will be made. Can be restricted to type. */


script_command *script_commands_create(struct script_session *s);
int command_compare(const void *a, const void *b);
void script_commands_print(FILE *f, const struct script_command *commands);
size_t script_commands_size(const script_command *commands);
void script_print_command_tree(FILE *f, const struct script_command *commands);
script_command_status script_execute_command(struct script_session *s, const char *cmd);
script_command_status script_execute_command_argv(struct script_session *s, const script_command *commands, int argc, char **argv);
void script_command_not_found(const char *cmd, const script_command *c);
const script_command *script_command_find(const script_command *commands, const char *cmd_string);
size_t script_command_print_possible_matches_if_ambiguous(const script_command *commands, const char *cmd_string);
script_command_status script_set_var(struct script_session *s, jibal_config_var *var, int, char * const *);
script_command_status script_show_var(struct script_session *s, jibal_config_var *var, int, char * const *);
script_command_status script_enable_var(struct script_session *s, jibal_config_var *var, int, char * const *);
script_command_status script_disable_var(struct script_session *s, jibal_config_var *var, int, char * const *);
script_command_status script_add_detector(struct script_session *s, int argc, char * const *argv);
script_command_status script_add_fit_range(struct script_session *s, int argc, char * const *argv);
script_command_status script_add_reaction(struct script_session *s, int argc, char * const *argv);
script_command_status script_add_reactions(struct script_session *s, int argc, char * const *argv);
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
script_command_status script_simulate(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_aperture(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_ion(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_detector(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_detector_aperture(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_detector_calibration(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_detector_foil(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_sample(struct script_session *s, int argc, char * const *argv);

#endif //JABS_SCRIPT_COMMAND_H
