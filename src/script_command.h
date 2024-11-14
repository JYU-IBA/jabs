/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2024 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */


#ifndef JABS_SCRIPT_COMMAND_H
#define JABS_SCRIPT_COMMAND_H

#include "script_generic.h"

#ifdef __cplusplus
extern "C" {
#endif

const char *script_command_status_to_string(script_command_status status);

script_command *script_command_new(const char *name, const char *help_text, int val, int argc_min, script_command_status (*f)(script_session *, int, char *const *)); /* Allocates new command that doesn't do anything. Can return var. */
int script_command_set_var(script_command *c, const jibal_config_var *var);
void script_command_free(script_command *c);
void script_commands_free(script_command *head);

script_command *script_command_list_find_tail(script_command *head);
script_command *script_command_list_merge(script_command *left, script_command *right);
script_command *script_command_list_merge_sort(script_command *head); /* Performs a merge sort of command list, sorts in alphabetical order by command name. Does not sort subcommands. */
script_command *script_command_list_append(script_command *head, script_command *cmd);
void script_command_list_add_command(script_command **head, script_command *c_new);
script_command *script_command_list_from_vars_array(const jibal_config_var *vars, jibal_config_var_type type); /* vars is an array, must be "null-terminated" (name of last element is NULL pointer). Deep copy will be made. Can be restricted to type. */

script_command *script_commands_create(script_session *s);
script_command *script_commands_sort_all(script_command *head); /* Sorts all commands (tree) so that each subcommand (branch) is sorted by script_command_list_merge_sort() */
void script_commands_print(const script_command *commands);
void script_print_command_tree(const script_command *commands);
script_command_status script_execute_command(script_session *s, const char *cmd);
script_command_status script_execute_command_argv(script_session *s, const script_command *commands, int argc, char **argv);
void script_command_not_found(const char *cmd, const script_command *c);
const script_command *script_command_find(const script_command *commands, const char *cmd_string);
size_t script_command_print_possible_matches_if_ambiguous(const script_command *commands, const char *cmd_string);
script_command_status script_set_var(script_session *s, jibal_config_var *var, int, char * const *);
script_command_status script_show_var(script_session *s, jibal_config_var *var, int, char * const *);
script_command_status script_set_charge(script_session *s, int argc, char *const *argv);
script_command_status script_set_detector_val(script_session *s, int val, int argc, char *const *argv);
script_command_status script_set_detector_calibration_val(script_session *s, int val, int argc, char *const *argv);
script_command_status script_set_fit_val(script_session *s, int val, int argc, char *const *argv);
script_command_status script_set_simulation_val(script_session *s, int val, int argc, char *const *argv);
script_command_status script_enable_var(script_session *s, jibal_config_var *var, int, char * const *);
script_command_status script_disable_var(script_session *s, jibal_config_var *var, int, char * const *);
script_command_status script_add_detector_default(script_session *s, int argc, char * const *argv);
script_command_status script_add_fit_range(script_session *s, int argc, char * const *argv);
script_command_status script_add_reaction(script_session *s, int argc, char * const *argv);
script_command_status script_add_reactions(script_session *s, int argc, char * const *argv);
script_command_status script_exit(script_session *s, int argc, char * const *argv);
script_command_status script_fit(script_session *s, int argc, char * const *argv);
script_command_status script_help(script_session *s, int argc, char * const *argv);
script_command_status script_help_version(script_session *s, int argc, char *const *argv);
script_command_status script_help_commands(script_session *s, int argc, char *const *argv);
#ifdef JABS_PLUGINS
script_command_status script_identify_plugin(script_session *s, int argc, char * const *argv);
script_command_status script_load_reaction_plugin(script_session *s, int argc, char *const *argv);
#endif
script_command_status script_load_experimental(script_session *s, int argc, char *const *argv);
script_command_status script_load_reference(script_session *s, int argc, char *const *argv);
script_command_status script_load_reaction(script_session *s, int argc, char *const *argv);
script_command_status script_load_roughness(script_session *s, int argc, char *const *argv);
script_command_status script_load_sample(script_session *s, int argc, char *const *argv);
script_command_status script_load_script(script_session *s, int argc, char *const *argv);
script_command_status script_remove_reaction(script_session *s, int argc, char * const *argv);
script_command_status script_reset(script_session *s, int argc, char * const *argv);
script_command_status script_reset_detectors(script_session *s, int argc, char * const *argv);
script_command_status script_reset_experimental(script_session *s, int argc, char * const *argv);
script_command_status script_reset_reference(script_session *s, int argc, char *const *argv);
script_command_status script_reset_fit(script_session *s, int argc, char * const *argv);
script_command_status script_reset_reactions(script_session *s, int argc, char * const *argv);
script_command_status script_reset_sample(script_session *s, int argc, char * const *argv);
script_command_status script_reset_stopping(script_session *s, int argc, char * const *argv);
script_command_status script_roi(script_session *s, int argc, char * const *argv);
script_command_status script_save_calibrations(script_session *s, int argc, char * const *argv);
script_command_status script_save_sample(script_session *s, int argc, char * const *argv);
script_command_status script_save_simulation(script_session *s, int argc, char * const *argv);
script_command_status script_save_spectra(script_session *s, int argc, char * const *argv);
script_command_status script_show_aperture(script_session *s, int argc, char * const *argv);
script_command_status script_show_calc_params(script_session *s, int argc, char * const *argv);
script_command_status script_show_detector(script_session *s, int argc, char * const *argv);
script_command_status script_show_fit(script_session *s, int argc, char * const *argv);
script_command_status script_show_fit_variables(script_session *s, int argc, char * const *argv);
script_command_status script_show_fit_ranges(script_session *s, int argc, char * const *argv);
script_command_status script_show_fit_correlation(script_session *s, int argc, char *const *argv);
script_command_status script_show_reactions(script_session *s, int argc, char * const *argv);
script_command_status script_show_sample(script_session *s, int argc, char * const *argv);
script_command_status script_show_sample_profile(script_session *s, int argc, char * const *argv);
script_command_status script_show_simulation(script_session *s, int argc, char * const *argv);
script_command_status script_show_stopping(script_session *s, int argc, char * const *argv);
script_command_status script_simulate(script_session *s, int argc, char * const *argv);
script_command_status script_set_channeling_yield(script_session *s, int argc, char *const *argv);
script_command_status script_set_channeling_slope(script_session *s, int argc, char *const *argv);
script_command_status script_set_aperture(script_session *s, int argc, char * const *argv);
script_command_status script_set_ion(script_session *s, int argc, char * const *argv);
script_command_status script_set_detector(script_session *s, int argc, char * const *argv);
script_command_status script_set_detector_aperture(script_session *s, int argc, char * const *argv);
script_command_status script_set_detector_calibration(script_session *s, int argc, char * const *argv);
script_command_status script_set_detector_foil(script_session *s, int argc, char * const *argv);
script_command_status script_set_detector_calibration_poly(script_session *s, int argc, char *const *argv);
script_command_status script_set_detector_name(script_session *s, int argc, char * const *argv);
script_command_status script_set_sample(script_session *s, int argc, char * const *argv);
script_command_status script_set_stopping(script_session *s, int argc, char * const *argv);
script_command_status script_set_verbosity(script_session *s, int argc, char * const *argv);
script_command_status script_test_reference(script_session *s, int argc, char * const *argv);
script_command_status script_test_roi(script_session *s, int argc, char * const *argv);
script_command_status script_split_sample_elements(script_session *s, int argc, char * const *argv);
script_command_status script_cwd(script_session *s, int argc, char * const *argv);
script_command_status script_cd(script_session *s, int argc, char * const *argv);
script_command_status script_idf2jbs(script_session *s, int argc, char * const *argv);
script_command_status script_kinematics(script_session *s, int argc, char * const *argv);
#ifdef __cplusplus
}
#endif
#endif //JABS_SCRIPT_COMMAND_H
