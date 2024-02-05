/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#include <stdlib.h>
#include <string.h>
#include "jabs_debug.h"
#include "plugin.h"
#include "generic.h"

jabs_plugin *jabs_plugin_open(const char *filename) {
    jabs_plugin_type (*typef)(void);
    const char *(*namef)(void);
    const char *(*versionf)(void);
    if(!filename) {
        return NULL;
    }
    void *handle = dlopen(filename, RTLD_LAZY);
    if(!handle) {
        return NULL;
    }
    jabs_plugin *plugin = malloc(sizeof(jabs_plugin));
    if(!plugin) {
        return NULL;
    }
    plugin->handle = handle;
    plugin->filename = strdup(filename);
    namef = (const char *(*)(void)) dlsym(handle, "name");
    versionf = (const char *(*)(void)) dlsym(handle, "version");
    typef = (jabs_plugin_type (*)(void)) dlsym(handle, "type");
    if(!namef || !versionf || !typef) {
        DEBUGSTR("Plugin does not have name or version symbol.");
        jabs_plugin_close(plugin);
        return NULL;
    }
    plugin->type = typef();
    plugin->name = strdup_non_null(namef());
    plugin->version = strdup_non_null(versionf());
    return plugin;
}

void jabs_plugin_close(jabs_plugin *plugin) {
    if(!plugin) {
        return;
    }
    dlclose(plugin->handle);
    free(plugin->filename);
    free(plugin->name);
    free(plugin->version);
    free(plugin);
}

const char *jabs_plugin_type_name(const jabs_plugin *plugin) {
    if(!plugin) {
        return NULL;
    }
    return jabs_plugin_type_string(plugin->type);
}

const char *jabs_plugin_type_string(jabs_plugin_type type) {
    switch(type) {
        case JABS_PLUGIN_NONE:
            return "none";
        case JABS_PLUGIN_CS:
            return "cs";
        case JABS_PLUGIN_SPECTRUM_READER:
            return "spectrumreader";
        default:
            return "unknown";
    }
}

jabs_plugin_reaction *jabs_plugin_reaction_init(const jabs_plugin *plugin, const jibal_isotope *isotopes, const jibal_isotope *incident, const jibal_isotope *target, int *argc, char * const **argv) {
    if(!plugin) {
        DEBUGSTR("Plugin is NULL. Can't init.");
        return NULL;
    }
    if(plugin->type != JABS_PLUGIN_CS) {
        return NULL;
    }
    jabs_plugin_reaction *(*initf)(const jibal_isotope *isotopes, const jibal_isotope *incident, const jibal_isotope *target, int *argc, char * const **argv);
    initf = (jabs_plugin_reaction *(*)(const jibal_isotope *, const jibal_isotope *, const jibal_isotope *, int *, char * const **)) dlsym(plugin->handle, "reaction_init");
    if(!initf) {
        DEBUGSTR("No reaction_init in plugin.\n");
        return NULL;
    }
    jabs_plugin_reaction *r = initf(isotopes, incident, target, argc, argv);
    return r;
}

void jabs_plugin_reaction_free(const jabs_plugin *plugin, jabs_plugin_reaction *reaction) {
    if(!plugin || !reaction) {
        return;
    }
    void *(*freef)(jabs_plugin_reaction *reaction);
    freef = (void *(*)(jabs_plugin_reaction *)) dlsym(plugin->handle, "reaction_free");
    if(!freef) {
        return;
    }
    freef(reaction);
}
