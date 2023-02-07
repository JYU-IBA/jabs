/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */


#ifndef JABS_SCRIPT_COMMAND_H
#define JABS_SCRIPT_COMMAND_H

#include "script_generic.h"

const char *script_command_status_to_string(script_command_status status);

script_command *script_command_new(const char *name, const char *help_text, int val, script_command_status (*f)(struct script_session *, int, char * const *)); /* Allocates new command that doesn't do anything. Can return var. */
int script_command_set_var(script_command *c, jibal_config_var_type type, void *variable, const jibal_option *option_list);
void script_command_free(script_command *c);
void script_commands_free(script_command *head);

script_command *script_command_list_find_tail(script_command *head);
script_command *script_command_list_merge(script_command *left, script_command *right);
script_command *script_command_list_merge_sort(script_command *head); /* Performs a merge sort of command list, sorts in alphabetical order by command name. Does not sort subcommands. */
script_command *script_command_list_append(script_command *head, script_command *cmd);
void script_command_list_add_command(script_command **head, script_command *c_new);
script_command *script_command_list_from_command_array(const script_command *commands); /* commands is an array, must be "null-terminated" (name of last element is NULL pointer). Deep copy will be made. */
script_command *script_command_list_from_vars_array(const jibal_config_var *vars, jibal_config_var_type type); /* vars is an array, must be "null-terminated" (name of last element is NULL pointer). Deep copy will be made. Can be restricted to type. */

script_command *script_commands_create(struct script_session *s);
script_command *script_commands_sort_all(script_command *head); /* Sorts all commands (tree) so that each subcommand (branch) is sorted by script_command_list_merge_sort() */
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
script_command_status script_set_detector_val(struct script_session *s, int val, int argc, char *const *argv);
script_command_status script_set_detector_calibration_val(struct script_session *s, int val, int argc, char *const *argv);
script_command_status script_set_fit_val(struct script_session *s, int val, int argc, char *const *argv);
script_command_status script_set_simulation_val(struct script_session *s, int val, int argc, char *const *argv);
script_command_status script_enable_var(struct script_session *s, jibal_config_var *var, int, char * const *);
script_command_status script_disable_var(struct script_session *s, jibal_config_var *var, int, char * const *);
script_command_status script_add_detector(struct script_session *s, int argc, char * const *argv);
script_command_status script_add_fit_range(struct script_session *s, int argc, char * const *argv);
script_command_status script_add_reaction(struct script_session *s, int argc, char * const *argv);
script_command_status script_add_reactions(struct script_session *s, int argc, char * const *argv);
script_command_status script_exit(struct script_session *s, int argc, char * const *argv);
script_command_status script_fit(struct script_session *s, int argc, char * const *argv);
script_command_status script_help(struct script_session *s, int argc, char * const *argv);
script_command_status script_help_version(struct script_session *s, int argc, char *const *argv);
script_command_status script_help_commands(script_session *s, int argc, char *const *argv);
#ifdef JABS_PLUGINS
script_command_status script_identify_plugin(struct script_session *s, int argc, char * const *argv);
script_command_status script_load_reaction_plugin(script_session *s, int argc, char *const *argv);
#endif
script_command_status script_load_experimental(struct script_session *s, int argc, char *const *argv);
script_command_status script_load_reference(struct script_session *s, int argc, char *const *argv);
script_command_status script_load_reaction(struct script_session *s, int argc, char *const *argv);
script_command_status script_load_roughness(script_session *s, int argc, char *const *argv);
script_command_status script_load_sample(struct script_session *s, int argc, char *const *argv);
script_command_status script_load_script(struct script_session *s, int argc, char *const *argv);
script_command_status script_remove_reaction(struct script_session *s, int argc, char * const *argv);
script_command_status script_reset(struct script_session *s, int argc, char * const *argv);
script_command_status script_reset_detectors(struct script_session *s, int argc, char * const *argv);
script_command_status script_reset_experimental(struct script_session *s, int argc, char * const *argv);
script_command_status script_reset_reference(script_session *s, int argc, char *const *argv);
script_command_status script_reset_fit(struct script_session *s, int argc, char * const *argv);
script_command_status script_reset_reactions(struct script_session *s, int argc, char * const *argv);
script_command_status script_reset_sample(struct script_session *s, int argc, char * const *argv);
script_command_status script_reset_stopping(struct script_session *s, int argc, char * const *argv);
script_command_status script_roi(struct script_session *s, int argc, char * const *argv);
script_command_status script_save_bricks(struct script_session *s, int argc, char * const *argv);
script_command_status script_save_calibrations(struct script_session *s, int argc, char * const *argv);
script_command_status script_save_sample(struct script_session *s, int argc, char * const *argv);
script_command_status script_save_spectra(struct script_session *s, int argc, char * const *argv);
script_command_status script_show_aperture(struct script_session *s, int argc, char * const *argv);
script_command_status script_show_calc_params(script_session *s, int argc, char * const *argv);
script_command_status script_show_detector(struct script_session *s, int argc, char * const *argv);
script_command_status script_show_fit(struct script_session *s, int argc, char * const *argv);
script_command_status script_show_fit(struct script_session *s, int argc, char * const *argv);
script_command_status script_show_fit_variables(struct script_session *s, int argc, char * const *argv);
script_command_status script_show_fit_ranges(struct script_session *s, int argc, char * const *argv);
script_command_status script_show_reactions(struct script_session *s, int argc, char * const *argv);
script_command_status script_show_sample(struct script_session *s, int argc, char * const *argv);
script_command_status script_show_simulation(struct script_session *s, int argc, char * const *argv);
script_command_status script_show_stopping(script_session *s, int argc, char * const *argv);
script_command_status script_simulate(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_channeling_yield(script_session *s, int argc, char *const *argv);
script_command_status script_set_channeling_slope(script_session *s, int argc, char *const *argv);
script_command_status script_set_aperture(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_ion(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_detector(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_detector_aperture(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_detector_calibration(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_detector_foil(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_detector_calibration_poly(struct script_session *s, int argc, char *const *argv);
script_command_status script_set_sample(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_stopping(struct script_session *s, int argc, char * const *argv);
script_command_status script_test_reference(struct script_session *s, int argc, char * const *argv);
script_command_status script_test_roi(struct script_session *s, int argc, char * const *argv);
script_command_status script_split_sample_elements(struct script_session *s, int argc, char * const *argv);
script_command_status script_cwd(struct script_session *s, int argc, char * const *argv);
script_command_status script_cd(struct script_session *s, int argc, char * const *argv);
script_command_status script_idf2jbs(struct script_session *s, int argc, char * const *argv);
#endif //JABS_SCRIPT_COMMAND_H
