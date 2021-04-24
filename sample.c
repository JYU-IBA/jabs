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
        if(x < sample->ranges[i + 1].x) {
            depth.x = x;
            depth.i = i;
            return depth;
        }
    }
    depth.i = sample->n_ranges-1;
    depth.x = sample->ranges[sample->n_ranges-1].x;
    return depth;
}

extern inline double depth_diff(depth a, depth b);

depth depth_add(const sample *sample, const depth in, double dx) {
    depth out;
    out.x = in.x + dx;
    size_t i = in.i;
    if(dx >= 0.0) {
        while(i < sample->n_ranges - 1) {
            if(out.x < sample->ranges[i+1].x) {
                //fprintf(stderr, "POS dx=%g tfu, %g tfu < %g tfu\n", dx/C_TFU, out.x/C_TFU, sample->cranges[i+1]/C_TFU);
                out.i = i;
                return out;
            }
            i++;
        }
        out.i = i;
#ifdef DEBUG
        if(i == sample->n_ranges-1) {
            fprintf(stderr, "Whoop! Time for the last depth bin! X = %g tfu, last known depth %g tfu.\n", out.x/C_TFU, sample->ranges[sample->n_ranges-1].x/C_TFU);
        }
#endif
        return out;
    } else { /* negative */
        while(i > 0) {
            if(out.x < sample->ranges[i-1].x) {
                fprintf(stderr, "NEG dx=%g tfu, %g tfu < %g tfu, i=%zu (from %zu, %g tfu)\n", dx/C_TFU, out.x/C_TFU, sample->ranges[i+1].x/C_TFU, i, in.i, in.x/C_TFU);
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
extern inline double *sample_conc_bin(const sample *s, size_t i_range, size_t i_isotope);
extern inline double *sample_model_conc_bin(const sample_model *sm, size_t i_range, size_t i_material);

extern inline int depth_is_almost_inside(double x, double low, double high) { /* Almost is good enough for us! */
    static const double abs_tol = DEPTH_TOLERANCE; /* we consider x to be in range [low, high] with this tolerance */
    return (x >= low-abs_tol && x <= high+abs_tol);
}

size_t get_range_bin(const sample *s, double x, size_t *range_hint) {
    int lo, mi, hi;
    if(range_hint) {
        if(*range_hint < s->n_ranges-1 && depth_is_almost_inside(x, s->ranges[*range_hint].x, s->ranges[*range_hint+1].x)) { /* TODO: add a bit of floating point "relative accuracy is enough" testing here */
            return *range_hint;
        } else { /* Shouldn't happen */
            fprintf(stderr, "WARNING!!! FALSE RANGE HINTING at depth = %g tfu. Hint was %lu (pointer %prob). ", x/C_TFU, *range_hint, (void *)range_hint);
            if(*range_hint >= s->n_ranges) {
                fprintf(stderr, "This is unacceptable because %lu should be < %lu.\n", *range_hint, s->n_ranges);
            } else {
                fprintf(stderr, "This is unacceptable because %g should be >= %g and <= %g.\n", x/C_TFU, s->ranges[*range_hint].x/C_TFU, s->ranges[*range_hint+1].x/C_TFU);
            }
        }
    } else { /* Shouldn't happen. Not wrong usage as such. */
        fprintf(stderr, "Mild warning: no range hinting, depth = %g tfu\n", x/C_TFU);
    }
    hi = s->n_ranges;
    lo = 0;
    while (hi - lo > 1) {
        mi = (hi + lo) / 2;
        if (x >= s->ranges[mi].x) {
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
    assert(x >= s->ranges[i_range].x);
    assert(x <= s->ranges[i_range+1].x);
    size_t i = i_range * s->n_isotopes + i_isotope;
    if(s->ranges[i_range+1].x - s->ranges[i_range].x == 0) /* Zero width. Return value of left side. */
        return s->cbins[i];
    return s->cbins[i] + ((s->cbins[i+s->n_isotopes] - s->cbins[i])/(s->ranges[i_range+1].x - s->ranges[i_range].x)) * (x - s->ranges[i_range].x);
}

int get_concs(const sample *s, const depth depth, double *out) {
    size_t i_range = depth.i;
    double x = depth.x;
    assert(i_range < s->n_ranges-1);
#if 0
    /* Due to floating point issues it is possible that depth.x is outside ranges[i] and ranges[i+1]. Note that depth.i should be respected even in this case. */
    assert(x >= s->ranges[i_range].x);
    assert(x <= s->ranges[i_range+1].x);
#endif
    double *bins_low = &s->cbins[i_range * s->n_isotopes];
    double *bins_high = &s->cbins[(i_range+1) * s->n_isotopes];
    if(s->ranges[i_range].x == s->ranges[i_range+1].x) {
        memcpy(out, bins_low, s->n_isotopes*sizeof(double));
        return 1;
    } else {
        double dx = (s->ranges[i_range+1].x - s->ranges[i_range].x);
        double deltax = (x - s->ranges[i_range].x);
        for (size_t i = 0; i < s->n_isotopes; i++) {
            out[i] = bins_low[i] + ((bins_high[i] - bins_low[i])/dx) * deltax;
        }
        return 2;
    }
}

sample *sample_alloc(size_t n_isotopes, size_t n_ranges) {
    sample *s = malloc(sizeof(sample));
    if(!s)
        return NULL;
    s->n_ranges = n_ranges;
    s->n_isotopes = n_isotopes;
    s->isotopes = calloc(n_isotopes, sizeof(jibal_isotope *));
    s->cbins = calloc( n_ranges * n_isotopes, sizeof(double));
    s->ranges = malloc(n_ranges*sizeof(struct sample_range));
    return s;
}

sample *sample_from_layers(jibal_layer * const *layers, size_t n_layers) {
    size_t i, j, k;
    sample *s = malloc(sizeof(sample));
    s->n_isotopes = 0;
    s->n_ranges = 2*n_layers;
    s->ranges = malloc(s->n_ranges*sizeof(struct sample_range));
    size_t i_isotope, n_isotopes=0;
    for(i = 0; i < n_layers; i++) { /* Turn layers to point-by-point profiles first by setting the X-axis (ranges) values and counting the number of isotopes */
        const jibal_layer *layer = layers[i];
#ifdef DEBUG
        fprintf(stderr, "Layer %lu/%lu. Thickness %g tfu\n", i+1, n_layers, layer->thickness/C_TFU);
        jibal_material_print(stderr, layer->material);
#endif
        s->ranges[2*i].x = i?s->ranges[2*i-1].x:0.0;
        s->ranges[2*i+1].x = s->ranges[2*i].x + (layer->thickness > 0.0 ? layer->thickness : 0.0); /* negative thickness... */
        s->ranges[2*i].rough.x = 0.0;
        s->ranges[2*i+1].rough.x = layer->roughness;
        for (j = 0; j < layer->material->n_elements; ++j) {
            n_isotopes += layer->material->elements[j].n_isotopes;
        }
    }
#ifdef DEBUG
    for (i = 0; i < s->n_ranges; i++) {
        fprintf(stderr, "ranges[%lu]  = %g\n", i, s->ranges[i].x);
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

sample_model *sample_model_alloc(size_t n_materials, size_t n_ranges) {
    sample_model *sm = malloc(sizeof(sample_model));
    if(!sm)
        return NULL;
    sm->type = SAMPLE_MODEL_NONE;
    sm->n_ranges = n_ranges;
    sm->n_materials = n_materials;
    sm->materials = calloc(n_materials, sizeof(jibal_material *));
    sm->cbins = calloc( n_ranges * n_materials, sizeof(double));
    sm->ranges = malloc(n_ranges*sizeof(struct sample_range));
    return sm;
}

sample_model *sample_model_copy(const sample_model *sm) {
    if(!sm)
        return NULL;
    sample_model *sm_copy = sample_model_alloc(sm->n_materials, sm->n_ranges);
    if(!sm_copy)
        return NULL;
//    memcpy(sm_copy->materials, )
    /* TODO */
    return sm_copy;
}

sample_model *sample_model_split_elements(const sample_model *sm) {
    sample_model *out = sample_model_alloc(sample_model_element_count(sm), sm->n_ranges);
    out->type = sm->type;
#ifdef DEBUG
    fprintf(stderr, "Splitting sample model with %zu ranges and %zu materials to one with %zu materials (=elements).\n", sm->n_ranges, sm->n_materials, out->n_materials);
#endif
    size_t i = 0; /* Material index in output */
    for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
        for(size_t i_elem = 0; i_elem < sm->materials[i_mat]->n_elements; i_elem++) {
            jibal_element *e = jibal_element_copy(&sm->materials[i_mat]->elements[i_elem], 0);
            out->materials[i] = malloc(sizeof(jibal_material));
            jibal_material *mat = out->materials[i];
            mat->elements = malloc(sizeof(jibal_element));
            mat->elements[0] = *e;
            mat->name = strdup(e->name);
            mat->n_elements = 1;
            mat->concs = malloc(sizeof(double));
            mat->concs[0] = 1.0;
            for(size_t i_range = 0; i_range < sm->n_ranges; i_range++) {
                *sample_model_conc_bin(out, i_range, i) += *sample_model_conc_bin(sm, i_range, i_mat) * sm->materials[i_mat]->concs[i_elem];
            }
            i++;
        }
    }
    for(size_t i_range = 0; i_range < sm->n_ranges; i_range++) {
        out->ranges[i_range] = sm->ranges[i_range];
    }
    return out;
}

size_t sample_model_element_count(const sample_model *sm) {
    size_t n_elements = 0;
    for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
        n_elements += sm->materials[i_mat]->n_elements;
    }
    return n_elements;
}

sample_model *sample_model_to_point_by_point(const sample_model *sm) { /* Converts a layered model to a point-by-point, allocates space and doesn't share data with the original (deep copy) */
    sample_model *sm_out = NULL;
    if(sm->type == SAMPLE_MODEL_LAYERED) { /* From layer model, make two points from one layer with the same concentration, i.e. no concentration gradient */
        sm_out = sample_model_alloc(sm->n_materials, sm->n_ranges*2);
        for(size_t i = 0; i < sm->n_materials; i++) {
            sm_out->materials[i] = jibal_material_copy(sm->materials[i]);
        }
        for(size_t i = sm->n_ranges; i--;) {
            sm_out->ranges[2*i] = sm->ranges[i];
            sm_out->ranges[2*i+1] = sm->ranges[i];
            for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
                *(sample_model_conc_bin(sm_out, 2*i, i_mat)) = *(sample_model_conc_bin(sm, i, i_mat));
                *(sample_model_conc_bin(sm_out, 2*i+1, i_mat)) = *(sample_model_conc_bin(sm, i, i_mat));
            }
        }
        for(size_t i = 0; i < sm_out->n_ranges; i += 2) {
            sm_out->ranges[i].rough.x = 0.0; /* First points are not rough, the second point carries this information */
            if(i) {
                sm_out->ranges[i].x = sm_out->ranges[i-1].x; /* In layer models, the first point has the same depth as the second point of previous layer */
                sm_out->ranges[i+1].x += sm_out->ranges[i].x; /* And the second point is naturally deeper by thickness, which was stored as "depth" temporarily */
            } else {
                sm_out->ranges[i].x = 0.0;
            }
        }
    }
    return sm_out;
}

sample *sample_from_sample_model(const sample_model *sm) {
    if(!sm)
        return NULL;
#ifdef DEBUG
    fprintf(stderr, "Sample model is type %i, it has %zu materials and %zu ranges.\n", sm->type, sm->n_materials, sm->n_ranges);
#endif
    sample_model *sm_copy;
    if(sm->type == SAMPLE_MODEL_LAYERED) {
        sm_copy = sample_model_to_point_by_point(sm);
        if(!sm_copy)
            return NULL;
        sm = sm_copy;
    } else {
        sm_copy = NULL;
    }
    size_t n_isotopes = 0;
    for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
        for(size_t i_elem = 0; i_elem < sm->materials[i_mat]->n_elements; i_elem++) {
            n_isotopes += sm->materials[i_mat]->elements[i_elem].n_isotopes;
        }
    }
    size_t i = 0;
    const jibal_isotope **isotopes = malloc(n_isotopes * sizeof(jibal_isotope *));
    for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
        for(size_t i_elem = 0; i_elem < sm->materials[i_mat]->n_elements; i_elem++) {
            for(size_t i_isotope = 0; i_isotope < sm->materials[i_mat]->elements[i_elem].n_isotopes; i_isotope++) {
                assert(i < n_isotopes);
                isotopes[i] = sm->materials[i_mat]->elements[i_elem].isotopes[i_isotope];
                i++;
            }
        }
    }
    sample *s = malloc(sizeof(sample));
    s->n_isotopes = n_isotopes;
    s->n_ranges = sm->n_ranges;
    s->isotopes = isotopes;
    sample_sort_and_remove_duplicate_isotopes(s);
    n_isotopes = s->n_isotopes;
    s->ranges = malloc(s->n_ranges * sizeof(struct sample_range));
    s->cbins = calloc(s->n_ranges * s->n_isotopes, sizeof(double));
    memcpy(s->ranges, sm->ranges, sizeof (struct sample_range) * sm->n_ranges);

    for(i = 0; i < s->n_isotopes; i++) { /* Build table of isotopic concentrations by looping over all isotopes in all elements in all materials and ranges*/
        for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
            for(size_t i_elem = 0; i_elem < sm->materials[i_mat]->n_elements; i_elem++) {
                for(size_t i_isotope = 0; i_isotope < sm->materials[i_mat]->elements[i_elem].n_isotopes; i_isotope++) {
                    if(s->isotopes[i] != sm->materials[i_mat]->elements[i_elem].isotopes[i_isotope])
                        continue;
                    for(size_t i_range = 0; i_range < sm->n_ranges; i_range++) {
                        *(sample_conc_bin(s, i_range, i)) += *(sample_model_conc_bin(sm, i_range, i_mat))  * sm->materials[i_mat]->elements[i_elem].concs[i_isotope] * sm->materials[i_mat]->concs[i_elem];
                    }
                }
            }
        }
    }
#ifdef DEBUG
    sample_print(stderr, s, 0);
#endif
    free(sm_copy);
    return s;
}

void sample_model_print(FILE *f, const sample_model *sm) {
    if(!sm)
        return;
    switch(sm->type) {
        case SAMPLE_MODEL_NONE:
            fprintf(stderr, "Sample model is none.\n");
            return;
            break;
        case SAMPLE_MODEL_POINT_BY_POINT:
            fprintf(f, "       depth");
            break;
        case SAMPLE_MODEL_LAYERED:
            fprintf(f, "       thick");
            break;
    }
    fprintf(f, "        rough");
    for (size_t i = 0; i < sm->n_materials; i++) {
        fprintf(f, " %8s", sm->materials[i]->name);
    }
    fprintf(f, "\n");
    for (size_t i = 0; i < sm->n_ranges; i++) {
        fprintf(f, "%12.3lf", sm->ranges[i].x/C_TFU);
        fprintf(f, " %12.3lf", sm->ranges[i].rough.x/C_TFU);
        for (size_t j = 0; j < sm->n_materials; j++) {
            fprintf(f, " %8.4lf", *sample_model_conc_bin(sm, i, j) * 100.0);
        }
        fprintf(f, "\n");
    }
}

sample_model *sample_model_from_file(jibal *jibal, const char *filename) {
    FILE *in;
    if(!filename)
        return NULL;
    in = fopen(filename, "r");
    if(!in)
        return NULL;
    sample_model *sm = malloc(sizeof(sample_model));
    sm->n_ranges = 0;
    sm->n_materials = 0;
    sm->ranges = NULL;
    sm->materials = NULL;
    sm->cbins = NULL;
    sm->type = SAMPLE_MODEL_NONE;

    char *line = NULL;
    size_t line_size = 0;
    size_t lineno = 0;

    size_t i_depth = 0, i_rough = 0;

    int headers = 1;

    while(getline(&line, &line_size, in) > 0) {
        lineno++;
        line[strcspn(line, "\r\n")] = 0; /* Strips all kinds of newlines! */
        if(strlen(line) >= 1 && *line == '#') /* Comment */
            continue;
        char *line_split = line;
        char *col;
        size_t n = 0; /* Number of columns on this row */
        size_t i_material = 0;
        while((col = strsep(&line_split, " \t")) != NULL) {
            if(*col == '\0') {/* Multiple separators are treated as one */
                continue;
            }
            if(headers) { /* Headers */
                if(strncmp(col, "rough", 5) == 0) {
                    i_rough = n;
                }
                else if(strncmp(col, "thick", 5) == 0) {
                    i_depth = n;
                    sm->type = SAMPLE_MODEL_LAYERED;
                }
                else if(strncmp(col, "depth", 5) == 0) {
                    i_depth = n;
                    sm->type = SAMPLE_MODEL_POINT_BY_POINT;
                } else {
                    sm->materials = realloc(sm->materials, sizeof (jibal_material *) * (sm->n_materials + 1)); /* TODO: check allocation */
                    sm->materials[sm->n_materials] = jibal_material_create(jibal->elements, col);
                    if(!sm->materials[sm->n_materials]) {
                        fprintf(stderr, "Could not create material %s, ignoring it!\n", col);
                    } else {
                        sm->n_materials++;
                    }
                }
                n++;
                continue;
            }
            if(n == 0) {
                sm->ranges = realloc(sm->ranges, sizeof(sample_range) * (sm->n_ranges + 1));
                sm->cbins = realloc(sm->cbins, sizeof(double) * (sm->n_ranges + 1) * sm->n_materials);
                assert(sm->cbins && sm->ranges);
                sm->ranges[sm->n_ranges].x = 0.0;
                sm->ranges[sm->n_ranges].rough.x = 0.0;
                sm->n_ranges++;
#ifdef DEBUG
                fprintf(stderr, "Sample from file: NEW POINT, n = %zu, n_ranges = %zu, n_materials=%zu\n", n, sm->n_ranges, sm->n_materials);
#endif
                i_material = 0;
                /* TODO: invalidate layer if not loaded properly */
            }
            sample_range *r = &sm->ranges[sm->n_ranges-1];

            double x = strtod(col, NULL);
            if(n == i_depth) {
                r->x += x*C_TFU; /* Will be corrected later in case of a layer model*/
                /* TODO: for p by p profile, check depth monotonicity */
            } else if (n == i_rough) {
                r->rough.x = x*C_TFU;
            } else if(i_material < sm->n_materials) {
                *(sample_model_conc_bin(sm, sm->n_ranges-1, i_material)) = x;
                i_material++;
            }
            n++;
        }
        headers = 0;
    }

    for(size_t i_range = 0; i_range < sm->n_ranges; i_range++) { /* Normalize concs */
        double sum = 0.0;
        for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
            sum += *(sample_model_conc_bin(sm, i_range, i_mat));
        }
        if(sum == 0.0)
            continue;
        for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
            *(sample_model_conc_bin(sm, i_range, i_mat)) /= sum;
        }
    }
    return sm;
}

void sample_model_free(sample_model *sm) {
    if(!sm)
        return;

    for(size_t i = 0; i < sm->n_materials; i++) {
        jibal_material_free(sm->materials[i]);
    }
    free(sm->ranges);
    free(sm->materials);
    free(sm->cbins);
    free(sm);
}

sample_model *sample_model_from_argv(jibal *jibal, int argc, char **argv) {
    sample_model *sm = malloc(sizeof(sample_model));
    sm->type = SAMPLE_MODEL_LAYERED;
    size_t n = 0;
    sm->n_ranges = 0;
    sm->n_materials = 0;
    sm->materials = NULL;
    sm->ranges = NULL;

    while (argc >= 2) {
        if(sm->n_ranges == n) {
            if(n == 0) {
                n = 8;
            } else {
                n *= 2;
            }
#ifdef DEBUG
            fprintf(stderr, "(Re)allocating space for up to %lu ranges.\n", n);
#endif
            sm->ranges = realloc(sm->ranges, n*sizeof(sample_range));
            sm->materials = realloc(sm->materials, n*sizeof(jibal_material *));
            if(!sm->ranges)
                return NULL;
        }
        if(sm->n_ranges && strcmp(argv[0], "rough") == 0) {
            sm->ranges[sm->n_ranges - 1].rough.x = jibal_get_val(jibal->units, UNIT_TYPE_LAYER_THICKNESS, argv[1]);
        } else {
            sm->materials[sm->n_ranges] = jibal_material_create(jibal->elements, argv[0]);
            sm->ranges[sm->n_ranges].x = jibal_get_val(jibal->units, UNIT_TYPE_LAYER_THICKNESS, argv[1]);
            sm->n_ranges++;
            sm->n_materials++;
        }
        argc -= 2;
        argv += 2;
    }
    sm->cbins = calloc(sm->n_ranges * sm->n_ranges, sizeof(double));
    for(size_t i = 0; i < sm->n_ranges; i++) {
        *sample_model_conc_bin(sm, i, i) = 1.0;
    }
    return sm;
}

sample *sample_copy(const sample *s_in) {
    sample *s_out = sample_alloc(s_in->n_isotopes, s_in->n_ranges);
    if(!s_out)
        return NULL;
    memcpy(s_out->isotopes, s_in->isotopes, sizeof(jibal_isotope *) * s_out->n_isotopes);
    memcpy(s_out->ranges, s_in->ranges, sizeof (struct sample_range) * s_out->n_ranges);
    memcpy(s_out->cbins, s_in->cbins, sizeof (double) * s_out->n_isotopes * s_out->n_ranges);
    return s_out;
}

void sample_areal_densities_print(FILE *f, const sample *sample, int print_isotopes) {
    fprintf(f, "AREAL D(tfu)             ");
    double sum = 0.0;
    for (size_t i = 0; i < sample->n_isotopes; i++) {
        for (size_t j = 1; j < sample->n_ranges; j++) {
            double thickness = (sample->ranges[j].x - sample->ranges[j-1].x);
            sum += 0.5 * ( *(sample_conc_bin(sample, j, i)) + *(sample_conc_bin(sample, j-1, i))) * thickness;
        }
        if (print_isotopes || i == sample->n_isotopes-1 || sample->isotopes[i]->Z != sample->isotopes[i+1]->Z) {
            fprintf(f, " %8.2lf", sum/C_TFU);
            sum = 0.0;
        }
    }
    fprintf(f, "\n");
}


void sample_print(FILE *f, const sample *sample, int print_isotopes) {
    fprintf(f, "  DEPTH(tfu)   ROUGH(tfu)");
    int Z = 0;
    for (size_t i = 0; i < sample->n_isotopes; i++) {
        if(print_isotopes) {
            fprintf(f, " %8s", sample->isotopes[i]->name);
        } else if(Z != sample->isotopes[i]->Z){
            const char *s = sample->isotopes[i]->name;
            while(*s >= '0' && *s <= '9') {s++;} /* Skip numbers, e.g. 28Si -> Si */
            fprintf(f, " %8s", s);
            Z = sample->isotopes[i]->Z; /* New element */
        }
    }
    fprintf(f, "\n");
    for (size_t i = 0; i < sample->n_ranges; i++) {
        fprintf(f, "%12.3lf", sample->ranges[i].x/C_TFU);
        fprintf(f, " %12.3lf", sample->ranges[i].rough.x/C_TFU);
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
    free(sample->ranges);
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
    return sample->ranges[i].x;
}

void sample_sort_isotopes(sample *sample) {
    qsort(sample->isotopes, sample->n_isotopes, sizeof(jibal_isotope *), &isotope_compar);
}

void sample_sort_and_remove_duplicate_isotopes(sample *s) {
    sample_sort_isotopes(s);
    size_t i_dst = 0;
    for(size_t i_src = 0; i_src < s->n_isotopes; i_src++) {
        if(i_src && s->isotopes[i_src] == s->isotopes[i_src-1])
            continue;
        s->isotopes[i_dst] = s->isotopes[i_src];
        i_dst++;
    }
    s->n_isotopes = i_dst;
    s->isotopes = realloc(s->isotopes, sizeof (jibal_isotope *) * s->n_isotopes); /* Frees unused memory */
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
