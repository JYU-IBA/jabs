/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2024 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_SCRIPT_GENERIC_H
#define JABS_SCRIPT_GENERIC_H

#include <time.h>
#include <jibal.h>
#include "defaults.h"
#include "fit.h"

#ifdef __cplusplus
extern "C" {
#endif

struct script_command; /* Declaration is required, because script_session has struct script_command pointer and script_command has script_session pointer */

typedef struct script_file {
    FILE *f;
    char *filename;
    char *line;
    size_t line_size;
    size_t lineno;
    int interactive;
} script_file;

typedef struct script_session {
    jibal *jibal;
    struct fit_data *fit;
    int (*fit_iter_callback)(fit_stats stats);
    double start, end; /* Time */
    script_file *files[SCRIPT_FILES_NESTED_MAX];
    size_t file_depth;
    struct script_command *commands;
    size_t i_det_active; /* Used by "set detector" */
    int Z_active; /* Used by "set detector calibration" */
} script_session;

#define SCRIPT_COMMAND_SUCCESS (0) /* Anything above and including zero is success */
#define SCRIPT_COMMAND_FAILURE (-1)
#define SCRIPT_COMMAND_NOT_FOUND (-2)
#define SCRIPT_COMMAND_EXIT (-3)
#define SCRIPT_COMMAND_EOF (-4)
#define SCRIPT_COMMAND_RESET (-5)

typedef int script_command_status; /* Script commands should return negative on error (see defines above) and number of arguments (zero or positive) consumed successfully */

typedef struct script_command {
    char *name;
    script_command_status (*f)(script_session *, int, char * const *); /* Function to process argument vectors. Called if it is non-NULL, even if subcommands exist! Return value is important.*/
    script_command_status (*f_var)(script_session *, jibal_config_var *var, int, char * const *); /* Function to process variables (from subcommands). */
    script_command_status (*f_val)(script_session *, int, int, char * const *); /* Function to process vals (from subcommands). */
    char *help_text; /* Short help text */
    struct script_command *subcommands;
    jibal_config_var *var;
    int val;
    struct script_command *next;
    int argc_min; /* Minimum number of arguments to run f. Set to zero to handle argc < argc_min in f. */
} script_command;

struct help_topic {
    const char *name;
    const char *help_text;
};

#ifdef __cplusplus
}
#endif
#endif // JABS_SCRIPT_GENERIC_H
