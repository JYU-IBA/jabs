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

extern inline double *sample_conc_bin(const sample *s, int i_range, int i_isotope);

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

double get_conc(sim_workspace *ws, const sample *s, double x, size_t i_isotope) {
    size_t i_range, i;

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
    size_t i_range, i;

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
    if(s->cranges[i_range] == s->cranges[i_range+1]) {
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

sample *sample_from_layers(jibal_layer * const *layers, size_t n_layers) {
    size_t i, j, k;
    sample *s = malloc(sizeof(sample));
    s->n_isotopes = 0;
    s->n_ranges = 2*n_layers;
    s->cranges = malloc(s->n_ranges*sizeof(double));
    size_t i_isotope, n_isotopes=0;
    for(i = 0; i < n_layers; i++) {
        const jibal_layer *layer = layers[i];
#ifdef DEBUG
        fprintf(stderr, "Layer %lu/%lu. Thickness %g tfu\n", i+1, n_layers, layer->thickness/C_TFU);
        jibal_material_print(stderr, layer->material);
#endif
        s->cranges[2*i] = i?s->cranges[2*i-1]:0.0;
        s->cranges[2*i+1] = s->cranges[2*i] + (layer->thickness > 0.0 ? layer->thickness : 0.0); /* negative thickness... */
        for (j = 0; j < layer->material->n_elements; ++j) {
            n_isotopes += layer->material->elements[j].n_isotopes;
        }
    }
#ifdef DEBUG
    for (i = 0; i < s->n_ranges; i++) {
        fprintf(stderr, "ranges[%lu]  = %g\n", i, s->cranges[i]);
    }
    fprintf(stderr, "Total %lu isotopes and %lu ranges\n", n_isotopes, s->n_ranges);
#endif
    s->isotopes = calloc(n_isotopes, sizeof(jibal_isotope *));
    i_isotope = 0;
    for(i = 0; i < n_layers; i++) { /* Set isotope pointers */
        const jibal_layer *layer = layers[i];
        for (j = 0; j < layer->material->n_elements; ++j) {
            for (k = 0; k < layer->material->elements[j].n_isotopes; k++) {
                assert(i_isotope < n_isotopes);
                s->isotopes[i_isotope] = layer->material->elements[j].isotopes[k];
                i_isotope++;
            }
        }
    }

#ifndef NO_LAYER_REMOVE_DUPLICATES

#ifdef DEBUG
    for(i = 0; i < n_isotopes; i++) {
        if(s->isotopes[i])
            fprintf(stderr, "%lu: %s\n", i, s->isotopes[i]->name);
    }
    fprintf(stderr, "Removing duplicates!\n");
#endif
    for(i = 0; i < n_isotopes; i++) { /* Set all duplicate isotope pointers to NULL */
        if(s->isotopes[i] == NULL)
            continue;
        for (j = i+1; j < n_isotopes; j++) {
            if(s->isotopes[i] == s->isotopes[j]) {
                s->isotopes[j] = NULL;
            }
        }
    }
#ifdef DEBUG
    for(i = 0; i < n_isotopes; i++) {
        if(s->isotopes[i]) {
            fprintf(stderr, "%lu: %s\n", i, s->isotopes[i]->name);
        }
    }
#endif
    for(i = 0; i < n_isotopes; i++) { /* Move NULL pointers to the end of array and calculate new size */
        if(s->isotopes[i] != NULL)
            continue;
#ifdef DEBUG
        fprintf(stderr, "shuffle i=%lu\n", i);
#endif
        for (j = i; j < n_isotopes; j++) {
            s->isotopes[j] = s->isotopes[j+1];
        }
        i--;
        n_isotopes--;
    }
#ifdef DEBUG
    fprintf(stderr, "%lu non-duplicate isotopes:\n", n_isotopes);
    for (i = 0; i < n_isotopes; i++) {
        fprintf(stderr, "%lu: %s\n", i, s->isotopes[i]->name);
    }
#endif
#endif // NO_LAYER_REMOVE_DUPLICATES


    s->cbins = calloc( s->n_ranges * n_isotopes, sizeof(double));
    s->n_isotopes = n_isotopes;
    /* TODO: it is possible to save a bit of memory by reallocing s->cbins and s->isotopes to match the new size (n_isotopes) */

    for(i_isotope = 0; i_isotope < n_isotopes; i_isotope++) {
        for (i = 0; i < n_layers; i++) {
            const jibal_layer *layer = layers[i];
            for (j = 0; j < layer->material->n_elements; j++) {
                if (layer->material->concs[j] < 0.0)
                    layer->material->concs[j] = 0.0001; /* TODO: fitting robustification */
            }
            jibal_material_normalize(layer->material);
            for (j = 0; j < layer->material->n_elements; j++) {
                jibal_element *element = &layer->material->elements[j];
                for (k = 0; k < element->n_isotopes; k++) {
                    if(layer->material->elements[j].isotopes[k] == s->isotopes[i_isotope]) {
                        s->cbins[(2 * i) * s->n_isotopes + i_isotope] = element->concs[k] * layer->material->concs[j];
                        s->cbins[(2 * i + 1) * s->n_isotopes + i_isotope] = element->concs[k] * layer->material->concs[j];
                    }
                }
            }
        }
    }
    return s;
}


void sample_print(FILE *f, const sample *sample) {
    fprintf(f, "DEPTH(tfu) ");
    for (size_t i = 0; i < sample->n_isotopes; i++) {
        fprintf(f, "%8s ", sample->isotopes[i]->name);
    }
    fprintf(f, "\n");
    for (size_t i = 0; i < sample->n_ranges; i++) {
        fprintf(f, "%10.3lf", sample->cranges[i]/C_TFU);
        for (size_t j = 0; j < sample->n_isotopes; j++) {
            fprintf(f, " %8.4lf", sample->cbins[i * sample->n_isotopes + j]*100.0);
        }
        fprintf(f, "\n");
    }
#ifdef PRINT_SAMPLE_MODEL
    fprintf(f, "\n");
    for (i = 0; i < 1000; ++i) {
        double d = 10.0*C_TFU*i;
        fprintf(f, "%10.3lf", d/C_TFU);
        for (j = 0; j < sample->n_isotopes; j++) {
            double x = get_conc(sample, d, j);
            fprintf(f, " %8.4lf", x*100.0);
        }
        fprintf(f, "\n");
    }
#endif
}

void sample_free(sample *sample) {
    if(!sample)
        return;
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
