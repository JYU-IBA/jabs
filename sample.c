/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */
#include <stdlib.h>
#include <assert.h>
#include "sample.h"


inline double *sample_conc_bin(const sample *s, int i_range, int i_isotope) {
    return s->cbins + i_range * s->n_isotopes + i_isotope;
}

int get_range_bin(const sample *s, double x) {
    int lo, mi, hi;
#ifdef RANGE_PEDANTIC
    if(x < s->cranges[0] || x >= s->cranges[s->n_ranges-1]) { /* Out of bounds concentration is zero. Maybe the execution doesn't go here if all goes as planned, so this could be changed to an assert. */
        return -1;
    }
#else
    if(x < s->cranges[0])
        return 0;
    if(x >= s->cranges[s->n_ranges-1]) { /* Out of bounds concentration is zero. Maybe the execution doesn't go here if all goes as planned, so this could be changed to an assert. */
        return -1;
    }
#endif
    hi = s->n_ranges;
    lo = 0;
    while (hi - lo > 1) {
        mi = (hi + lo) / 2;
        if (x >= s->cranges[mi]) {
            lo = mi;
        } else {
            hi = mi;
        }
    }
    return lo;
}

double get_conc(sim_workspace *ws, const sample *s, double x, int i_isotope) {
    int i_range, i;

    if(x >= s->cranges[ws->i_range_accel] && x < s->cranges[ws->i_range_accel+1]) { /* Use saved bin. We'll probably receive a lot of repeated calls to the same bin, so this avoid cost of other range checking and binary searches. */
        i = ws->i_range_accel * s->n_isotopes + i_isotope;
        if(s->cranges[ws->i_range_accel] == s->cranges[ws->i_range_accel+1]) /* Constant */
            return s->cbins[i];
        else /* Linear interpolation */
            return s->cbins[i] + ((s->cbins[i+s->n_isotopes] - s->cbins[i])/(s->cranges[ws->i_range_accel+1] - s->cranges[ws->i_range_accel])) * (x - s->cranges[ws->i_range_accel]);
    }
    i_range = get_range_bin(s, x);
    ws->i_range_accel = i_range;
    if(i_range < 0) {
#ifdef DEBUG
        fprintf(stderr, "No depth range found for x = %.5lf\n", x/C_TFU);
#endif
        return 0.0;
    }
    assert(i_isotope < s->n_isotopes);
    i = i_range * s->n_isotopes + i_isotope;
    if(s->cranges[i_range+1] - s->cranges[i_range] == 0) /* Zero width. Return value of left side. */
        return s->cbins[i];
    return s->cbins[i] + ((s->cbins[i+s->n_isotopes] - s->cbins[i])/(s->cranges[i_range+1] - s->cranges[i_range])) * (x - s->cranges[i_range]);
}

int get_concs(sim_workspace *ws, const sample *s, double x, double *out) {
    int i_range, i;

    if(x >= s->cranges[ws->i_range_accel] && x < s->cranges[ws->i_range_accel+1]) { /* Use saved bin. */
        i_range = ws->i_range_accel;
    } else {
        i_range = get_range_bin(s, x);
        ws->i_range_accel = i_range;
    }
    if(i_range < 0) {
#ifdef DEBUG
        fprintf(stderr, "No depth range found for x = %.5lf\n", x/C_TFU);
#endif
        return 0;
    }
    double *bins_low = &s->cbins[i_range * s->n_isotopes];
    double *bins_high = &s->cbins[(i_range+1) * s->n_isotopes];
    if(s->cranges[i_range+1] - s->cranges[i_range] == 0) {
        for(i = 0; i < s->n_isotopes; i++) {
            out[i] = bins_low[i]; /* TODO: memcpy? */
        }
        return 1;
    } else {
        double dx = (s->cranges[i_range+1] - s->cranges[i_range]);
        double deltax = (x - s->cranges[i_range]);
        for (i = 0; i < s->n_isotopes; i++) {
            out[i] = bins_low[i] + ((bins_high[i] - bins_low[i])/dx) * deltax;
        }
        return 2;
    }
}

sample *sample_from_layers(jibal_layer **layers, int n_layers) {
    int i, j, k;
    sample *s = malloc(sizeof(sample));
    s->n_isotopes = 0;
    s->n_ranges = 2*n_layers;
    s->cranges = malloc(s->n_ranges*sizeof(double));
    int i_isotope, n_isotopes=0;
    for(i = 0; i < n_layers; i++) {
        jibal_layer *layer = layers[i];
#ifdef DEBUG
        fprintf(stderr, "Layer %i/%i. Thickness %g tfu\n", i+1, n_layers, layer->thickness/C_TFU);
        jibal_material_print(stderr, layer->material);
#endif
        s->cranges[2*i] = i?s->cranges[2*i-1]:0.0;
        s->cranges[2*i+1] = s->cranges[2*i] + layer->thickness;
        for (j = 0; j < layer->material->n_elements; ++j) {
            n_isotopes += layer->material->elements[j].n_isotopes;
        }
    }
#ifdef DEBUG
    for (i = 0; i < s->n_ranges; i++) {
        fprintf(stderr, "ranges[%i]  = %g\n", i, s->cranges[i]);
    }
    fprintf(stderr, "Total %i isotopes and %i ranges\n", n_isotopes, s->n_ranges);
#endif
    s->isotopes = calloc(n_isotopes, sizeof(jibal_isotope *));
    s->cbins = calloc( s->n_ranges * n_isotopes, sizeof(double));
    i_isotope = 0;
    s->n_isotopes = n_isotopes;
    for (i = 0; i < n_layers; i++) {
        jibal_layer *layer = layers[i];
        for (j = 0; j < layer->material->n_elements; ++j) {
            jibal_element *element = &layer->material->elements[j];
            for (k = 0; k < element->n_isotopes; k++) {
                //assert(i_isotope < n_isotopes);
                s->isotopes[i_isotope] = element->isotopes[k];
                s->cbins[(2 * i)*s->n_isotopes + i_isotope] = element->concs[k] * layer->material->concs[j];
                s->cbins[(2 * i + 1)*s->n_isotopes + i_isotope] = element->concs[k] * layer->material->concs[j];
                i_isotope++;
            }
        }
    }

    return s;
}


void sample_print(FILE *f, const sample *sample) {
    int i,j;
    fprintf(stderr, "DEPTH(tfu) ");
    for (i = 0; i < sample->n_isotopes; i++) {
        fprintf(stderr, "%8s ", sample->isotopes[i]->name);
    }
    fprintf(stderr, "\n");
    for (i = 0; i < sample->n_ranges; i++) {
        fprintf(stderr, "%10.3lf", sample->cranges[i]/C_TFU);
        for (j = 0; j < sample->n_isotopes; j++) {
            fprintf(stderr, " %8.4lf", sample->cbins[i * sample->n_isotopes + j]*100.0);
        }
        fprintf(stderr, "\n");
    }
#ifdef PRINT_SAMPLE_MODEL
    fprintf(stderr, "\n");
    for (i = 0; i < 1000; ++i) {
        double d = 10.0*C_TFU*i;
        fprintf(stderr, "%10.3lf", d/C_TFU);
        for (j = 0; j < sample->n_isotopes; j++) {
            double x = get_conc(sample, d, j);
            fprintf(stderr, " %8.4lf", x*100.0);
        }
        fprintf(stderr, "\n");
    }
#endif
}

void sample_free(sample *sample) {
    free(sample->cranges);
    free(sample->isotopes);
    free(sample->cbins);
    free(sample);
}

double sample_isotope_max_depth(const sample *sample, int i_isotope) {
    int i;
    for(i = sample->n_ranges - 1; i >= 0; i--) {
        double *c = sample_conc_bin(sample, i, i_isotope);
        if(*c > 0.0) /* TODO: other cutoff? */
            break;
    }
    return sample->cranges[i];
}
