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
#include <string.h>
#include "defaults.h"
#include "sample.h"

depth depth_seek(const sample *sample, double x) {
    depth depth;
    for(size_t i = 0; i < sample->n_ranges - 1; i++) {
        if(x < sample->cranges[i + 1]) {
            depth.x = x;
            depth.i = i;
            return depth;
        }
    }
    depth.i = sample->n_ranges-1;
    depth.x = sample->cranges[sample->n_ranges-1];
    return depth;
}

extern inline double depth_diff(depth a, depth b);

depth depth_add(const sample *sample, const depth in, double dx) {
    depth out;
    out.x = in.x + dx;
    size_t i = in.i;
    if(dx >= 0.0) {
        while(i < sample->n_ranges - 1) {
            if(out.x < sample->cranges[i+1]) {
                //fprintf(stderr, "POS dx=%g tfu, %g tfu < %g tfu\n", dx/C_TFU, out.x/C_TFU, sample->cranges[i+1]/C_TFU);
                out.i = i;
                return out;
            }
            i++;
        }
        out.i = i;
#ifdef DEBUG
        if(i == sample->n_ranges-1) {
            fprintf(stderr, "Whoop! Time for the last depth bin! X = %g tfu, last known depth %g tfu.\n", out.x/C_TFU, sample->cranges[sample->n_ranges-1]/C_TFU);
        }
#endif
        return out;
    } else { /* negative */
        while(i > 0) {
            if(out.x < sample->cranges[i-1]) {
                fprintf(stderr, "NEG dx=%g tfu, %g tfu < %g tfu, i=%zu (from %zu, %g tfu)\n", dx/C_TFU, out.x/C_TFU, sample->cranges[i+1]/C_TFU, i, in.i, in.x/C_TFU);
                out.i = i;
                return out;
            }
            i--;
        }
        out.i = i;
        return out;
    }
    /* No execution path reaches here */
}
extern inline double *sample_conc_bin(const sample *s, int i_range, int i_isotope);

extern inline int depth_is_almost_inside(double x, double low, double high) { /* Almost is good enough for us! */
    static const double abs_tol = DEPTH_TOLERANCE; /* we consider x to be in range [low, high] with this tolerance */
    return (x >= low-abs_tol && x <= high+abs_tol);
}

size_t get_range_bin(const sample *s, double x, size_t *range_hint) {
    int lo, mi, hi;
    if(range_hint) {
        if(*range_hint < s->n_ranges-1 && depth_is_almost_inside(x, s->cranges[*range_hint], s->cranges[*range_hint+1])) { /* TODO: add a bit of floating point "relative accuracy is enough" testing here */
            return *range_hint;
        } else { /* Shouldn't happen */
            fprintf(stderr, "WARNING!!! FALSE RANGE HINTING at depth = %g tfu. Hint was %lu (pointer %p). ", x/C_TFU, *range_hint, (void *)range_hint);
            if(*range_hint >= s->n_ranges) {
                fprintf(stderr, "This is unacceptable because %lu should be < %lu.\n", *range_hint, s->n_ranges);
            } else {
                fprintf(stderr, "This is unacceptable because %g should be >= %g and <= %g.\n", x/C_TFU, s->cranges[*range_hint]/C_TFU, s->cranges[*range_hint+1]/C_TFU);
            }
        }
    } else { /* Shouldn't happen. Not wrong usage as such. */
        fprintf(stderr, "Mild warning: no range hinting, depth = %g tfu\n", x/C_TFU);
    }
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

double get_conc(const sample *s, const depth depth, size_t i_isotope) {
    assert(i_isotope < s->n_isotopes);
    size_t i_range = depth.i;
    double x = depth.x;
    assert(i_range < s->n_ranges-1);
    assert(x >= s->cranges[i_range]);
    assert(x <= s->cranges[i_range+1]);
    size_t i = i_range * s->n_isotopes + i_isotope;
    if(s->cranges[i_range+1] - s->cranges[i_range] == 0) /* Zero width. Return value of left side. */
        return s->cbins[i];
    return s->cbins[i] + ((s->cbins[i+s->n_isotopes] - s->cbins[i])/(s->cranges[i_range+1] - s->cranges[i_range])) * (x - s->cranges[i_range]);
}

int get_concs(const sample *s, const depth depth, double *out) {
    size_t i_range = depth.i;
    double x = depth.x;
    assert(i_range < s->n_ranges-1);
    assert(x >= s->cranges[i_range]);
    assert(x <= s->cranges[i_range+1]);
    double *bins_low = &s->cbins[i_range * s->n_isotopes];
    double *bins_high = &s->cbins[(i_range+1) * s->n_isotopes];
    if(s->cranges[i_range] == s->cranges[i_range+1]) {
        memcpy(out, bins_low, s->n_isotopes*sizeof(double));
        return 1;
    } else {
        double dx = (s->cranges[i_range+1] - s->cranges[i_range]);
        double deltax = (x - s->cranges[i_range]);
        for (size_t i = 0; i < s->n_isotopes; i++) {
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
    for(i = 0; i < n_layers; i++) { /* Turn layers to point-by-point profiles first by setting the X-axis (ranges) values and counting the number of isotopes */
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

    sample_sort_isotopes(s); /* TODO: removing duplicates is easier after this step */

    for(i_isotope = 0; i_isotope < n_isotopes; i_isotope++) { /* Now the Y-values (bins) of the point-by-point profile */
        for (i = 0; i < n_layers; i++) {
            const jibal_layer *layer = layers[i];
            for (j = 0; j < layer->material->n_elements; j++) {
                if (layer->material->concs[j] < 0.0)
                    layer->material->concs[j] = DEPTH_TOLERANCE*10.0; /* TODO: fitting robustification */
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


void sample_print(FILE *f, const sample *sample, int print_isotopes) {
    fprintf(f, "DEPTH(tfu) ");
    int Z = 0;
    for (size_t i = 0; i < sample->n_isotopes; i++) {
        if(print_isotopes) {
            fprintf(f, "%8s ", sample->isotopes[i]->name);
        } else if(Z != sample->isotopes[i]->Z){
            const char *s = sample->isotopes[i]->name;
            while(*s >= '0' && *s <= '9') {s++;} /* Skip numbers, e.g. 28Si -> Si */
            fprintf(f, "%8s", s);
            Z = sample->isotopes[i]->Z; /* New element */
        }
    }
    fprintf(f, "\n");
    for (size_t i = 0; i < sample->n_ranges; i++) {
        fprintf(f, "%10.3lf", sample->cranges[i]/C_TFU);
        double sum = 0.0;
        for (size_t j = 0; j < sample->n_isotopes; j++) {
            sum += sample->cbins[i * sample->n_isotopes + j];
            if (print_isotopes || j == sample->n_isotopes-1 || sample->isotopes[j]->Z != sample->isotopes[j+1]->Z) { /* Last isotope or next isotope belongs to another element, print. */
                fprintf(f, " %8.4lf", sum * 100.0);
                sum = 0.0;
            }
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

void sample_sort_isotopes(sample *sample) {
    qsort(sample->isotopes, sample->n_isotopes, sizeof(jibal_isotope *), &isotope_compar);
}

int isotope_compar(const void *a, const void *b) {
    const jibal_isotope *isotope_a = *((const jibal_isotope **)a);
    const jibal_isotope *isotope_b = *((const jibal_isotope **)b);
    if(isotope_a->Z == isotope_b->Z) { /* Same Z, compare by A */
        return isotope_a->A - isotope_b->A;
    } else {
        return isotope_a->Z - isotope_b->Z;
    }
}
