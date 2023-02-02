/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

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
#include <jibal.h>
#include "roughness.h"

typedef struct sample_range {
    double x;
    roughness rough;
    double yield; /* Ad-hoc correction to scattering yield (cross section) */
    double bragg; /* Ad-hoc correction to stopping */
    double stragg; /* Ad-hoc correction to straggling */
    double density; /* Density of range (preferably of a layer) , used only to convert at./cm2 into nm if > 0.0 */
} sample_range;

typedef struct depth {
    double x;
    size_t i; /* index in sample */
} depth;

typedef enum sample_model_type {
    SAMPLE_MODEL_NONE = 0,
    SAMPLE_MODEL_POINT_BY_POINT = 1,
    SAMPLE_MODEL_LAYERED = 2
} sample_model_type;

typedef struct sample {
    int no_conc_gradients; /* TRUE (non zero) for samples based on layer models */
    size_t n_isotopes;
    size_t n_ranges;
    const jibal_isotope **isotopes; /* table, size is n_isotopes */
    sample_range *ranges; /* size is n_ranges */
    double *cbins; /* 2D-table: size is n_isotopes * n_ranges  */
    double thickness;
} sample;

typedef struct sample_model {
    sample_model_type type;
    size_t n_materials;
    size_t n_ranges;
    jibal_material **materials;
    sample_range *ranges;
    double *cbins; /* 2D-table: size is n_materials * n_ranges  */
} sample_model;

inline double *sample_model_conc_bin(const sample_model *sm, size_t i_range, size_t i_material) {
    return sm->cbins + i_range * sm->n_materials + i_material;
}
inline double *sample_conc_bin(const sample *s, size_t i_range, size_t i_isotope) {
    return s->cbins + i_range * s->n_isotopes + i_isotope;
}
sample_model *sample_model_alloc(size_t n_materials, size_t n_ranges);
int sample_model_sanity_check(const sample_model *sm);
void sample_model_renormalize(sample_model *sm);
void sample_renormalize(sample *sample);
sample_model *sample_model_split_elements(const struct sample_model *sm);
sample_model *sample_model_from_file(const jibal *jibal, const char *filename);
sample_model *sample_model_from_argv(const jibal *jibal, int *argc, char * const **argv);
sample_model *sample_model_from_string(const jibal *jibal, const char *str);
char *sample_model_to_string(const sample_model *sm);
sample_model *sample_model_to_point_by_point(const sample_model *sm);
size_t sample_model_element_count(const sample_model *sm);
void sample_model_free(sample_model *sm);
sample *sample_from_sample_model(const sample_model *sm);
int sample_model_print(const char *filename, const sample_model *sm);
size_t sample_model_number_of_rough_ranges(const sample_model *sm);
size_t sample_model_number_of_ranges_with_non_unity_corrections(const sample_model *sm);
size_t sample_model_number_of_range_with_non_zero_density(const sample_model *sm);
void sample_thickness_recalculate(sample *sample);

depth depth_seek(const sample *sample, double x);
inline double depth_diff(const depth a, const depth b) {
    return b.x - a.x;
}
double get_conc_interpolate(const sample *s, depth depth, size_t i_isotope);
inline double get_conc(const sample *s, const depth depth, size_t i_isotope) {
    return s->no_conc_gradients ? *sample_conc_bin(s, depth.i, i_isotope) : get_conc_interpolate(s, depth, i_isotope);
}
sample *sample_alloc(size_t n_isotopes, size_t n_ranges);
sample *sample_copy(const sample *sample); /* Deep copy */
double sample_mass_density_range(const sample *sample, size_t i_range);
double sample_thickness_in_nm_range(const sample *sample, size_t i_range);
double sample_areal_density_isotope_range(const sample *sample,  size_t i_isotope, size_t i_range);
void sample_areal_densities_print(FILE *f, const sample *sample, int print_isotopes);
int sample_print_thicknesses(const char *filename, const sample *sample);
int sample_print(const char *filename, const sample *sample, int print_isotopes);  /* If print_isotopes is non-zero print print_isotopes individually. Isotopes must be sorted by Z, e.g. with sample_sort_isotopes() */
void sample_free(sample *sample);
double sample_isotope_max_depth(const sample *sample, size_t i_isotope);
void sample_sort_isotopes(sample *sample); /* Sort isotopes. Note that you can screw up concentrations tables if you call this at the wrong time! */
void sample_sort_and_remove_duplicate_isotopes(sample *s);
int isotope_compar(const void *, const void *); /* compares jibal_isotopes */
void sample_range_copy(sample_range *dst, const sample_range *src);
#endif /* JABS_SAMPLE_H */
