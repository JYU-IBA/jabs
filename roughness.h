/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#ifndef JABS_ROUGHNESS_H
#define JABS_ROUGHNESS_H

typedef struct {
    double x;
    double prob;
} thick_prob;

typedef struct {
    thick_prob *p;
    size_t n;
} thick_prob_dist;

typedef enum roughness_model {
    ROUGHNESS_NONE = 0,
    ROUGHNESS_GAUSSIAN = 1, /* TODO: implement! */
    ROUGHNESS_GAMMA = 2
} roughness_model;

typedef struct roughness {
    roughness_model model;
    double x; /* Amount, model specific? */
    size_t n; /* Number of spectra (e.g. for gamma roughness) */
} roughness;

thick_prob_dist *thickness_probability_table_gen(double thickness, double sigma, size_t n);
void thickness_probability_table_free(thick_prob_dist *tpd);
double thickness_gamma_pdf(double x, double thickness, double sigma);
#endif // JABS_ROUGHNESS_H
