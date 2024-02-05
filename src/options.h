/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_OPTIONS_H
#define JABS_OPTIONS_H

#include <jibal.h>

typedef struct {
    int verbose;
    int interactive;
} cmdline_options;


#include "simulation.h"

void read_options(cmdline_options *cmd_opt, int *argc, char *const **argv);
cmdline_options *cmdline_options_init(void);
void cmdline_options_free(cmdline_options *cmd_opt);
const char *jabs_version(void); /* Returns "git describe" given version, fallback to jabs_version_simple() */
const char *jabs_version_simple(void); /* Returns directly the CMake project version */
void usage(void);
void greeting(int interactive);

#endif // JABS_OPTIONS_H
