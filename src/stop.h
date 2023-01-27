/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_STOP_H
#define JABS_STOP_H
#include "jibal_gsto.h"
#include "ion.h"
#include "sample.h"

typedef struct jabs_stop {
    jibal_gsto *gsto;
    gsto_stopping_type type;
    int nuclear_stopping_accurate;
    int rk4;
} jabs_stop;

depth stop_next_crossing(const ion *incident, const sample *sample, const depth *d_from);
depth stop_step(const jabs_stop *stop, const jabs_stop *stragg, ion *incident, const sample *sample, depth depth_before, double step);
double stop_sample(const jabs_stop *stop, const ion *incident, const sample *sample, depth depth, double E);
#endif // JABS_STOP_H
