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

typedef struct sample {
    size_t n_isotopes;
    size_t n_ranges;
    const jibal_isotope **isotopes; /* table, size is n_isotopes */
    double *cranges;
    double *crange_roughness;
    double *cbins; /* 2D-table: size is n_isotopes * n_ranges  */
} sample;

typedef struct depth {
    double x;
    size_t i; /* index in sample */
} depth;

depth depth_seek(const sample *sample, double x);
inline double depth_diff(const depth a, const depth b) {
    return b.x - a.x;
}
depth depth_add(const sample *sample, depth in, double dx); /* Create a new depth by adding (or substracting, if negative) thickness to depth "in") */
double get_conc(const sample *s, depth depth, size_t i_isotope);
int get_concs(const sample *s, depth depth, double *out);
size_t get_range_bin(const sample *s, double x, size_t *range_hint);

sample *sample_alloc(size_t n_isotopes, size_t n_ranges);
sample *sample_from_layers(jibal_layer * const *layers, size_t n_layers);
sample *sample_copy(const sample *sample); /* Deep copy */
void sample_areal_densities_print(FILE *f, const sample *sample, int print_isotopes);
void sample_print(FILE *f, const sample *sample, int print_isotopes); /* If print_isotopes is non-zero print print_isotopes individually. Isotopes must be sorted by Z, e.g. with sample_sort_isotopes() */
void sample_free(sample *sample);
double sample_isotope_max_depth(const sample *sample, int i_isotope);
inline double *sample_conc_bin(const sample *s, size_t i_range, size_t i_isotope) {
    return s->cbins + i_range * s->n_isotopes + i_isotope;
}

void sample_sort_isotopes(sample *sample); /* Sort isotopes. Note that you can screw up concentrations tables if you call this at the wrong time! */
int isotope_compar(const void *, const void *); /* compares jibal_isotopes */

#endif /* JABS_SAMPLE_H */
