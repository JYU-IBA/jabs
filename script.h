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
script_command_status script_show(script_session *s, int argc, char * const *argv);
script_command_status script_set(script_session *s, int argc, char * const *argv);
script_command_status script_add(script_session *s, int argc, char * const *argv);
script_command_status script_reset(script_session *s, int argc, char * const *argv);
script_command_status script_simulate(script_session *s, int argc, char * const *argv);
script_command_status script_fit(script_session *s, int argc, char * const *argv);
script_command_status script_save_spectra(script_session *s, int argc, char * const *argv);
script_command_status script_save_detector(script_session *s, int argc, char * const *argv);
script_command_status script_save_sample(script_session *s, int argc, char * const *argv);
script_command_status script_save(script_session *s, int argc, char * const *argv);
script_command_status script_remove(script_session *s, int argc, char * const *argv);
script_command_status script_roi(script_session *s, int argc, char * const *argv);
void script_command_not_found(const char *cmd, const script_command *parent);
int script_process(script_session *s, const char *filename);
int script_prepare_sim_or_fit(script_session *s);
int script_finish_sim_or_fit(script_session *s);
int script_get_detector_number(const simulation *sim, int *argc, char *const **argv, size_t *i_det);
static const struct script_command script_save_commands[] = {
        {"spectra",  &script_save_spectra,  "Save spectra.",  NULL},
        {"detector", &script_save_detector, "Save detector.", NULL},
        {"sample",   &script_save_sample,   "Save sample.",   NULL},
        {NULL, NULL, NULL,                                    NULL}
};

static const struct script_command script_commands[] = {
        {"help",     &script_help,     "Print help.",                                 NULL},
        {"show",     &script_show,     "Show information on things.",                 NULL},
        {"set",      &script_set,      "Set variables.",                              NULL},
        {"add",      &script_add,      "Add things.",                                 NULL},
        {"simulate", &script_simulate, "Run a simulation.",                           NULL},
        {"load",     &script_load,     "Load something.",                             NULL},
        {"reset",    &script_reset,    "Reset something.",                            NULL},
        {"fit",      &script_fit,      "Do a fit.",                                   NULL},
        {"roi",      &script_roi,      "Show information from a region of interest.", NULL},
        {"save",     NULL,     "Save something.",                             script_save_commands},
        {"remove",   &script_remove,   "Remove something",                            NULL},
        {"exit", NULL,                 "Exit.",                                       NULL},
        {"quit", NULL, NULL,                                                          NULL},
        {NULL,   NULL, NULL,                                                          NULL}
}; /* TODO: more commands... */
#endif // JABS_SCRIPT_H
