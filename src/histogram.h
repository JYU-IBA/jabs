/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See LICENSE.txt for the full license.

 */
#ifndef JABS_HISTOGRAM_H
#define JABS_HISTOGRAM_H

typedef struct {
    size_t n;
    double *range;
    double *bin;
} jabs_histogram; /* This is the same as gsl_histogram. */

/* These functions are re-implemented to avoid GSL error checking. */
jabs_histogram *jabs_histogram_clone(const jabs_histogram *h_orig); /* Doesn't do range checking */
void jabs_histogram_free(jabs_histogram *h);
jabs_histogram *jabs_histogram_alloc(size_t n);
void jabs_histogram_reset(jabs_histogram *h);
void jabs_histogram_scale(jabs_histogram *h, double scale);
inline double jabs_histogram_get(const jabs_histogram *h, size_t i) {return h->bin[i];}

/* These functions are JaBS additions */
int jabs_histogram_compare(const jabs_histogram *h1, const jabs_histogram *h2, size_t low, size_t high, double *out);
double jabs_histogram_roi(const jabs_histogram *h, size_t low, size_t high); /* low and high are both inclusive channel numbers*/
size_t jabs_histogram_channels_in_range(const jabs_histogram *h, size_t low, size_t high);
#endif // JABS_HISTOGRAM_H
