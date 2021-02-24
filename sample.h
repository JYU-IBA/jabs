/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_SAMPLE_H
#define JABS_SAMPLE_H

#include <stdio.h>
#include <jibal_masses.h>
#include <jibal_material.h>
#include <jibal_layer.h>


typedef struct {
    int n_isotopes;
    int n_ranges;
    const jibal_isotope **isotopes; /* table, size is n_isotopes */
    double *cranges;
    double *cbins; /* 2D-table: size is n_isotopes * n_ranges  */
} sample;

#include "simulation.h"

double get_conc(sim_workspace *ws, const sample *s, double x, int i);
int get_concs(sim_workspace *ws, const sample *s, double x, double *out);
int get_range_bin(const sample *s, double x);

sample *sample_from_layers(jibal_layer * const *layers, int n_layers);
void sample_print(FILE *f, const sample *sample);
void sample_free(sample *sample);
double sample_isotope_max_depth(const sample *sample, int i_isotope);
inline double *sample_conc_bin(const sample *s, int i_range, int i_isotope) {
    return s->cbins + i_range * s->n_isotopes + i_isotope;
}

#endif /* JABS_SAMPLE_H */
