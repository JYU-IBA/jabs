/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_TESTPLUGIN_H
#define JABS_TESTPLUGIN_H
#include "../src/plugin_interface.h"

const char *name(void);
const char *version(void);
jabs_plugin_type type(void);
jabs_plugin_reaction *reaction_init(const jibal_isotope *isotopes, const jibal_isotope *incident, const jibal_isotope *target, int *argc, char * const **argv);
void reaction_free(jabs_plugin_reaction *r);
double testplugin_cs(const struct jabs_plugin_reaction *r, double theta, double E);
#endif // JABS_TESTPLUGIN_H
