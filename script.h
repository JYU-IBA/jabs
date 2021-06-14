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


struct script_command {
    const char *name;
    int (*f)(struct fit_data *, jibal_config_var *, int, char * const *);
    const char *help_text; /* Short help text */
};

struct help_topic {
    const char *name;
    const char *help_text;
};

typedef struct script_session {
    jibal *jibal;
    struct fit_data *fit;
    jibal_config_var *vars;
} script_session;

script_session *script_session_init(jibal *jibal);
void script_session_free(script_session *s);

void script_print_commands(FILE *f);
int script_load(struct fit_data *fit, jibal_config_var *vars, int argc, char * const *argv);
int script_help(struct fit_data *fit, jibal_config_var *vars, int argc, char * const *argv);
int script_show(struct fit_data *fit, jibal_config_var *vars, int argc, char * const *argv);
int script_set(struct fit_data *fit, jibal_config_var *vars, int argc, char * const *argv);
int script_add(struct fit_data *fit, jibal_config_var *vars, int argc, char * const *argv);
int script_reset(struct fit_data *fit, jibal_config_var *vars, int argc, char * const *argv);
jibal_config_var *script_make_vars(struct fit_data *fit);
int script_simulate(struct fit_data *fit, jibal_config_var *vars, int argc, char * const *argv);
int script_fit(struct fit_data *fit_data, jibal_config_var *vars, int argc, char * const *argv);
int script_save(struct fit_data *fit_data, jibal_config_var *vars, int argc, char * const *argv);
int script_process(script_session *s, FILE *f);
int script_prepare_sim_or_fit(struct fit_data *fit_data);
#endif // JABS_SCRIPT_H
