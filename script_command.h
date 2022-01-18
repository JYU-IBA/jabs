/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021-2022 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

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

void script_print_commands(FILE *f, const struct script_command *commands);
void script_print_command_tree(FILE *f, const struct script_command *commands);
script_command_status script_execute_command(script_session *s, const char *cmd);
script_command_status script_execute_command_argv(script_session *s, const script_command *commands, int argc, char **argv);
void script_command_not_found(const char *cmd, const script_command *c);
const script_command *script_command_find(const script_command *commands, const char *cmd_string);
script_command_status script_set_boolean(script_session *s, const char *variable, int value);

script_command_status script_add_detector(script_session *s, int argc, char * const *argv);
script_command_status script_add_fit_range(script_session *s, int argc, char * const *argv);
script_command_status script_add_reaction(script_session *s, int argc, char * const *argv);
script_command_status script_add_reactions(script_session *s, int argc, char * const *argv);
script_command_status script_disable(script_session *s, int argc, char * const *argv);
script_command_status script_enable(script_session *s, int argc, char * const *argv);
script_command_status script_exit(script_session *s, int argc, char * const *argv);
script_command_status script_fit(script_session *s, int argc, char * const *argv);
script_command_status script_help(script_session *s, int argc, char * const *argv);
script_command_status script_load_detector(script_session *s, int argc, char *const *argv);
script_command_status script_load_experimental(script_session *s, int argc, char *const *argv);
script_command_status script_load_reaction(script_session *s, int argc, char *const *argv);
script_command_status script_load_sample(script_session *s, int argc, char *const *argv);
script_command_status script_load_script(script_session *s, int argc, char *const *argv);
script_command_status script_remove_reaction(script_session *s, int argc, char * const *argv);
script_command_status script_reset(script_session *s, int argc, char * const *argv);
script_command_status script_reset_detectors(script_session *s, int argc, char * const *argv);
script_command_status script_reset_experimental(script_session *s, int argc, char * const *argv);
script_command_status script_reset_fit_ranges(script_session *s, int argc, char * const *argv);
script_command_status script_reset_reactions(script_session *s, int argc, char * const *argv);
script_command_status script_reset_sample(script_session *s, int argc, char * const *argv);
script_command_status script_roi(script_session *s, int argc, char * const *argv);
script_command_status script_save_bricks(script_session *s, int argc, char * const *argv);
script_command_status script_save_detector(script_session *s, int argc, char * const *argv);
script_command_status script_save_sample(script_session *s, int argc, char * const *argv);
script_command_status script_save_spectra(script_session *s, int argc, char * const *argv);
script_command_status script_show_detector(script_session *s, int argc, char * const *argv);
script_command_status script_show_fit(script_session *s, int argc, char * const *argv);
script_command_status script_show_reactions(script_session *s, int argc, char * const *argv);
script_command_status script_show_sample(script_session *s, int argc, char * const *argv);
script_command_status script_show_simulation(script_session *s, int argc, char * const *argv);
script_command_status script_show_variables(script_session *s, int argc, char * const *argv);
script_command_status script_simulate(script_session *s, int argc, char * const *argv);
script_command_status script_set(script_session *s, int argc, char * const *argv);
script_command_status script_set_aperture(script_session *s, int argc, char * const *argv);
script_command_status script_set_beam(script_session *s, int argc, char * const *argv);
script_command_status script_set_ion(script_session *s, int argc, char * const *argv);
script_command_status script_set_detector(script_session *s, int argc, char * const *argv);
script_command_status script_set_sample(script_session *s, int argc, char * const *argv);
script_command_status script_set_variable(script_session *s, int argc, char * const *argv);

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
