/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_DES_H
#define JABS_DES_H
#include "sample.h"
#include "ion.h"
#include "stop.h"
#include "simulation.h"

typedef struct {
    depth d; /* Depth */
    double E; /* Energy */
    double S; /* Straggling */
} des; /* DES = Depth, Energy, Straggling */

typedef struct {
    des *t; /* array, n elements. E decreasing, d either increases (ion going deeper) or decreases. S can do whatever S does. */
    int depth_increases; /* TRUE if d increases (going deeper) */
    size_t n;
    size_t n_ranges;
    size_t *depth_interval_index; /* table, size same as number of sample ranges (as given in des_table_compute()). Array stores location i of t[i].d.x == sample->range[i_range].x  */
} des_table; /* Depth, energy, straggling table */

int des_table_realloc(des_table *dt, size_t n);
void des_table_free(des_table *dt);
des_table *des_table_compute(const jabs_stop *stop, const jabs_stop *stragg, const sim_calc_params *scp, const sample *sample, const ion *incident, depth depth_start, double emin);
size_t des_table_size(const des_table *dt);
des *des_table_element(const des_table *dt, size_t i);
void des_table_rebuild_index(des_table *dt); /* called by des_table_compute() after setting values to table and before any other function can be used */
void des_table_print(FILE *f, const des_table *dt);
depth des_table_find_depth(const des_table *dt, size_t *i_des, depth depth_prev, ion *incident); /* Returns depth at given incident->E, or the next layer boundary, starting search from i_des in DES table. Updates i_des, incident->E and ->S. */
inline const des *des_table_min_energy_bin(const des_table *dt) {return &(dt->t[dt->n - 1]);}
void des_set_ion(const des *des, ion *ion);
#endif // JABS_DES_H
