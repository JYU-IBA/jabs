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
} script_command;

struct help_topic {
    const char *name;
    const char *help_text;
};

typedef struct script_command_option { /* Command options are parsed by name, store an int and optionally a pointer to script command (c) or configuration variable (var). */
    const char *name;
    int val;
    const script_command *c;
    jibal_config_var *var;
} script_command_option; /* Note that name, c and var are all const pointers. No deep copies are made! */

script_command_option *options_from_commands(const script_command *commands);
script_command_option *options_from_jibal_config(jibal_config_var *vars);
script_command_option *options_append(script_command_option *opt_to, const script_command_option *opt_from);
void options_sort(script_command_option *opt);
int option_compare(const void *a, const void *b);
size_t options_size(const script_command_option *opt);
void options_print(FILE *f, const script_command_option *options);
char *options_list_matches(const script_command_option *options, const char *str);
int script_getopt(struct script_session *s, const script_command_option *options, int *argc, char *const **argv, script_command_status *status_out); /* Parses argument vector, finds script_command_option "c" and calls "c->f()" or sets c->var. */

void script_print_commands(FILE *f, const struct script_command *commands);
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
script_command_status script_set(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_aperture(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_beam(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_ion(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_detector(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_sample(struct script_session *s, int argc, char * const *argv);
script_command_status script_set_variable(struct script_session *s, int argc, char * const *argv);


static const struct script_command script_set_commands[] = {
        {"aperture", &script_set_aperture, "Set aperture.",               NULL},
        {"beam",     &script_set_beam,     "Set beam properties.",        NULL},
        {"detector", &script_set_detector, "Set detector properties.",    NULL},
        {"ion",      &script_set_ion,      "Set incident ion (isotope).", NULL},
        {"sample",   &script_set_sample,   "Set sample.",                 NULL},
        {"variable", &script_set_variable, "Set a variable",              NULL},
        {NULL, NULL, NULL,                                                NULL}
};

static const struct script_command script_load_commands[] = {
        {"detector",     &script_load_detector,     "Load (replace) a detector.",     NULL},
        {"experimental", &script_load_experimental, "Load an experimental spectrum.", NULL},
        {"script",       &script_load_script,       "Load (run) a script.",           NULL},
        {"sample",       &script_load_sample,       "Load a sample.",                 NULL},
        {"reaction",     &script_load_reaction,     "Load a reaction from R33 file.", NULL},
        {NULL, NULL, NULL,                                                            NULL}
};


static const struct script_command script_save_commands[] = {
        {"bricks",   &script_save_bricks,   "Save bricks.",   NULL},
        {"detector", &script_save_detector, "Save detector.", NULL},
        {"sample",   &script_save_sample,   "Save sample.",   NULL},
        {"spectra",  &script_save_spectra,  "Save spectra.",  NULL},
        {NULL, NULL, NULL,                                    NULL}
};

static const struct script_command script_add_commands[] = {
        {"detector",  &script_add_detector,  "Add a detector.",               NULL},
        {"fit_range", &script_add_fit_range, "Add a fit range",               NULL},
        {"reaction",  &script_add_reaction,  "Add a reaction.",               NULL},
        {"reactions", &script_add_reactions, "Add reactions (of some type).", NULL},
        {NULL, NULL, NULL,                                                    NULL}
};

static const struct script_command script_show_commands[] = {
        {"detector",   &script_show_detector,   "Show detector.",   NULL},
        {"fit",        &script_show_fit,        "Show fit.",        NULL},
        {"reactions",  &script_show_reactions,  "Show reactions,",  NULL},
        {"sample",     &script_show_sample,     "Show sample",      NULL},
        {"simulation", &script_show_simulation, "Show simulation.", NULL},
        {"variables",  &script_show_variables,  "Show variables.",  NULL},
        {NULL, NULL, NULL,                                          NULL}
};

static const struct script_command script_remove_commands[] = {
        {"reaction", &script_remove_reaction, "Remove reaction.", NULL},
        {NULL, NULL, NULL,                                        NULL}
};

static const struct script_command script_reset_commands[] = {
        {"detectors",    &script_reset_detectors,    "Reset detectors.",            NULL},
        {"experimental", &script_reset_experimental, "Reset experimental spectra.", NULL},
        {"fit_ranges",   &script_reset_fit_ranges,   "Reset fit ranges.",           NULL},
        {"reactions",    &script_reset_reactions,    "Reset reactions.",            NULL},
        {"sample",       &script_reset_sample,       "Reset sample.",               NULL},
        {NULL, NULL, NULL,                                                          NULL}
};

static const struct script_command script_commands[] = {
        {"add",    NULL,               "Add things.",                 script_add_commands},
        {"disable",  &script_disable,  "Set boolean variable to false.",              NULL},
        {"enable",   &script_enable,   "Set boolean variable to true.",               NULL},
        {"exit",     &script_exit,     "Exit.",                                       NULL},
        {"fit",      &script_fit,      "Do a fit.",                                   NULL},
        {"save",   NULL,               "Save something.",             script_save_commands},
        {"set",      &script_set,      "Set variables.",              script_set_commands},
        {"show",   NULL,               "Show information on things.", script_show_commands},
        {"help",     &script_help,     "Print help.",                                 NULL},
        {"load",   NULL,               "Load something.",             script_load_commands},
        {"remove", NULL,               "Remove something",            script_remove_commands},
        {"reset",    &script_reset,    "Reset something.",            script_reset_commands},
        {"roi",      &script_roi,      "Show information from a region of interest.", NULL},
        {"simulate", &script_simulate, "Run a simulation.",                           NULL},
        {NULL,     NULL, NULL,                                                        NULL},
};

#endif //JABS_SCRIPT_COMMAND_H
