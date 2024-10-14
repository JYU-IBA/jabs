/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2024 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#ifndef JABS_ROUGHNESS_H
#define JABS_ROUGHNESS_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double x;
    double prob;
} thick_prob;

typedef struct {
    thick_prob *p;
    size_t n;
    size_t modulo; /* cumulative product of number of thickness probability distributions. When we have two rough layers with n1 and n2, we need to simulate n1*n2 subspectra. */
    size_t i_range; /* depth bin */
} thick_prob_dist;

typedef enum roughness_model {
    ROUGHNESS_NONE = 0,
    ROUGHNESS_GAUSSIAN = 1, /* TODO: implement! */
    ROUGHNESS_GAMMA = 2,
    ROUGHNESS_FILE = 3
} roughness_model;

typedef struct roughness_file {
    char *filename;
    thick_prob_dist *tpd;
} roughness_file;

typedef struct roughness {
    roughness_model model;
    double x; /* Amount, model specific. */
    size_t n; /* Number of spectra to simulate (e.g. for gamma roughness), not used for roughness files. */
    roughness_file *file;
} roughness;


thick_prob_dist *thickness_probability_table_gamma(double thickness, double sigma, size_t n);
thick_prob_dist *thickness_probability_table_new(size_t n);
thick_prob_dist *thickness_probability_table_from_file(const char *filename);
void thickness_probability_table_normalize(thick_prob_dist *tpd);
double thickness_probability_table_areal_density(thick_prob_dist *tpd);
void thickness_probability_table_free(thick_prob_dist *tpd);
void thickness_probability_table_print(FILE *f, const thick_prob_dist *tpd);
thick_prob_dist *thickness_probability_table_realloc(thick_prob_dist *tpd, size_t n);
thick_prob_dist *thickness_probability_table_copy(const thick_prob_dist *tpd);
int roughness_reset(roughness *r);
int roughness_reset_if_below_tolerance(roughness *r);
int roughness_set_from_file(roughness *r, const char *filename);
roughness_file *roughness_file_read(const char *filename);
roughness_file *roughness_file_copy(const roughness_file *rf);
void roughness_file_free(roughness_file *rf);
double thickness_gamma_pdf(double x, double thickness, double sigma);
#ifdef __cplusplus
}
#endif
#endif // JABS_ROUGHNESS_H
