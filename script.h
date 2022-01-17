/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#ifndef JABS_SCRIPT_H
#define JABS_SCRIPT_H
#include <stdio.h>
#include "jabs.h"

typedef struct script_session {
    jibal *jibal;
    struct fit_data *fit;
    jibal_config_file *cf; /* a pseudo-configuration "file" for our session */
    char *output_filename; /* File name for automatic spectra saving */ /* TODO: multidetector! */
    char *bricks_out_filename; /* File name for automatic bricks saving */
    char *sample_out_filename; /* File name for automatic sample saving */
    char *detector_out_filename; /* File name for automatic detector saving */
    clock_t start, end;
} script_session;

typedef enum script_command_status {
    SCRIPT_COMMAND_SUCCESS = EXIT_SUCCESS,
    SCRIPT_COMMAND_FAILURE = EXIT_FAILURE,
    SCRIPT_COMMAND_NOT_FOUND = 101,
    SCRIPT_COMMAND_AMBIGUOUS = 102,
    SCRIPT_COMMAND_EXIT = 103
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

script_session *script_session_init(jibal *jibal, simulation *sim); /* sim can be NULL or a previously initialized sim can be given. Note that it will be free'd by script_session_free()! */
int script_session_reset_vars(script_session *s);
void script_session_free(script_session *s);

void script_print_commands(FILE *f, const struct script_command *commands);
jibal_config_var *script_make_vars(script_session *s);
script_command_status script_load(script_session *s, int argc, char * const *argv);
script_command_status script_help(script_session *s, int argc, char * const *argv);
script_command_status script_show_sample(script_session *s, int argc, char * const *argv);
script_command_status script_show_simulation(script_session *s, int argc, char * const *argv);
script_command_status script_show_fit(script_session *s, int argc, char * const *argv);
script_command_status script_show_detector(script_session *s, int argc, char * const *argv);
script_command_status script_show_reactions(script_session *s, int argc, char * const *argv);
script_command_status script_show_variables(script_session *s, int argc, char * const *argv);
script_command_status script_set(script_session *s, int argc, char * const *argv);
script_command_status script_add_reaction(script_session *s, int argc, char * const *argv);
script_command_status script_add_reactions(script_session *s, int argc, char * const *argv);
script_command_status script_add_detector(script_session *s, int argc, char * const *argv);
script_command_status script_add_fit_range(script_session *s, int argc, char * const *argv);
script_command_status script_reset(script_session *s, int argc, char * const *argv);
script_command_status script_simulate(script_session *s, int argc, char * const *argv);
script_command_status script_fit(script_session *s, int argc, char * const *argv);
script_command_status script_save_spectra(script_session *s, int argc, char * const *argv);
script_command_status script_save_detector(script_session *s, int argc, char * const *argv);
script_command_status script_save_sample(script_session *s, int argc, char * const *argv);

script_command_status script_remove(script_session *s, int argc, char * const *argv);
script_command_status script_roi(script_session *s, int argc, char * const *argv);
script_command_status script_exit(script_session *s, int argc, char * const *argv);
const script_command *script_command_find(const script_command *commands, const char *cmd_string); /* Returns pointer to command (if it is unambiguous) */
void script_command_not_found(const char *cmd, const script_command *parent);
int script_process(script_session *s, const char *filename);
int script_prepare_sim_or_fit(script_session *s);
int script_finish_sim_or_fit(script_session *s);
int script_get_detector_number(const simulation *sim, int *argc, char *const **argv, size_t *i_det);
script_command_status script_load_script(script_session *s, int argc, char *const *argv);
script_command_status script_load_sample(script_session *s, int argc, char *const *argv);
script_command_status script_load_detector(script_session *s, int argc, char *const *argv);
script_command_status script_load_experimental(script_session *s, int argc, char *const *argv);
script_command_status script_load_reaction(script_session *s, int argc, char *const *argv);
script_command_status script_remove_reaction(script_session *s, int argc, char * const *argv);

static const struct script_command script_load_commands[] = {
        {"detector",     &script_load_detector,     "Load (replace) a detector.",     NULL},
        {"experimental", &script_load_experimental, "Load an experimental spectrum.", NULL},
        {"script",       &script_load_script,       "Load (run) a script.",           NULL},
        {"sample",       &script_load_sample,       "Load a sample.",                 NULL},
        {"reaction",     &script_load_reaction,     "Load a reaction from R33 file.", NULL},
        {NULL, NULL, NULL,                                                            NULL}
};


static const struct script_command script_save_commands[] = {
        {"detector", &script_save_detector, "Save detector.", NULL},
        {"sample",   &script_save_sample,   "Save sample.",   NULL},
        {"spectra",  &script_save_spectra,  "Save spectra.",  NULL},
        {NULL, NULL, NULL, NULL}
};

static const struct script_command script_add_commands[] = {
        {"detector",   &script_add_detector,   "Add a detector.",   NULL},
        {"fit_range", &script_add_fit_range, "Add a fit range",               NULL},
        {"reaction",  &script_add_reaction,  "Add a reaction.",               NULL},
        {"reactions", &script_add_reactions, "Add reactions (of some type).", NULL},
        {NULL, NULL, NULL, NULL}
};

static const struct script_command script_show_commands[] = {
        {"detector",   &script_show_detector,   "Show detector.",   NULL},
        {"fit",        &script_show_fit,        "Show fit.",        NULL},
        {"reactions",  &script_show_reactions,  "Show reactions,",  NULL},
        {"sample",     &script_show_sample,     "Show sample",      NULL},
        {"simulation", &script_show_simulation, "Show simulation.", NULL},
        {"variables",  &script_show_variables,  "Show variables.",  NULL},
        {NULL, NULL, NULL,NULL}
};

static const struct script_command script_remove_commands[] = {
        {"reaction",   &script_remove_reaction,   "Remove reaction.",   NULL},
        {NULL, NULL, NULL,NULL}
};


static const struct script_command script_commands[] = {
        {"add",    NULL,               "Add things.",                 script_add_commands},
        {"exit",     &script_exit,     "Exit.",                                       NULL},
        {"fit",      &script_fit,      "Do a fit.",                                   NULL},
        {"save",   NULL,               "Save something.",             script_save_commands},
        {"set",      &script_set,      "Set variables.",                              NULL},
        {"show",   NULL,               "Show information on things.", script_show_commands},
        {"help",     &script_help,     "Print help.",                                 NULL},
        {"load",   NULL,               "Load something.",             script_load_commands},
        {"remove", NULL,               "Remove something",            script_remove_commands},
        {"reset",    &script_reset,    "Reset something.",                            NULL},
        {"roi",      &script_roi,      "Show information from a region of interest.", NULL},
        {"simulate", &script_simulate, "Run a simulation.",                           NULL},
        {NULL,     NULL, NULL,                                                        NULL},
};
#endif // JABS_SCRIPT_H
