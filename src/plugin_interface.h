/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_PLUGIN_INTERFACE_H
#define JABS_PLUGIN_INTERFACE_H
#include <jibal_masses.h>

typedef enum jabs_plugin_type {
    JABS_PLUGIN_NONE = 0,
    JABS_PLUGIN_CS = 1, /* Cross section evaluation plugin */
    JABS_PLUGIN_SPECTRUM_READER = 2 /* Spectrum reader plugin */
} jabs_plugin_type;

typedef struct jabs_plugin_reaction {
    const jibal_isotope *incident;
    const jibal_isotope *target;
    const jibal_isotope *product;
    const jibal_isotope *product_heavy;
    double (*cs)(const struct jabs_plugin_reaction *r, double theta, double E);
    double E_min, E_max;
    void *plugin_data; /* Pointer for the plugin to store reaction specific data in. */
} jabs_plugin_reaction;
#endif //JABS_PLUGIN_INTERFACE_H
