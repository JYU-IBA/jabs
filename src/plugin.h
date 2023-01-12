/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_PLUGIN_H
#define JABS_PLUGIN_H

#include <dlfcn.h> /* TODO: Windows? */
#include "plugin_interface.h"

typedef struct jabs_plugin {
    void *handle;
    char *filename;
    char *name; /* Symbol "name" from library */
    char *version; /* Symbol "version" from library */
    jabs_plugin_type type;
} jabs_plugin;

jabs_plugin *jabs_plugin_open(const char *filename);
void jabs_plugin_close(jabs_plugin *plugin);

const char *jabs_plugin_type_name(const jabs_plugin *plugin);
const char *jabs_plugin_type_string(jabs_plugin_type type);

jabs_plugin_reaction *jabs_plugin_reaction_init(const jabs_plugin *plugin, const jibal_isotope *incident, const jibal_isotope *target, int *argc, char * const **argv);
void jabs_plugin_reaction_free(const jabs_plugin *plugin, jabs_plugin_reaction *reaction);
#endif //JABS_PLUGIN_H
