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
    int (*f)(struct fit_data *, int, char * const *);
    const char *help_text; /* Short help text */
};

struct help_topic {
    const char *name;
    const char *help_text;
};

void script_print_commands(FILE *f);
int script_help(struct fit_data *fit, int argc, char * const *argv);
int script_show(struct fit_data *fit, int argc, char * const *argv);
int script_reset(struct fit_data *fit, int argc, char * const *argv);
int script_process(jibal *jibal, FILE *f);
#endif // JABS_SCRIPT_H
