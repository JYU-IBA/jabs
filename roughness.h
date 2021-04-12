/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

typedef struct {
    double x;
    double prob;
} thick_prob;

typedef struct {
    thick_prob *p;
    size_t n;
} thick_prob_dist;


thick_prob_dist *thickness_probability_table_gen(double thickness, double sigma, size_t n);
void thickness_probability_table_free(thick_prob_dist *tpd);
double thickness_gamma_pdf(double x, double thickness, double sigma);
