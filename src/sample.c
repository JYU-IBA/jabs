/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2022 Jaakko Julin

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
#include "generic.h"
#include "sample.h"
#include "message.h"
#include "win_compat.h"
#include "generic.h"

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
extern inline double *sample_conc_bin(const sample *s, size_t i_range, size_t i_isotope);
extern inline double *sample_model_conc_bin(const sample_model *sm, size_t i_range, size_t i_material);

double get_conc(const sample *s, const depth depth, size_t i_isotope) {
    assert(i_isotope < s->n_isotopes);
    size_t i_range = depth.i;
    double x = depth.x;
    assert(i_range < s->n_ranges-1);
#ifdef RANGE_PEDANTIC
    assert(x >= s->ranges[i_range].x);
    assert(x <= s->ranges[i_range+1].x);
#endif
    size_t i = i_range * s->n_isotopes + i_isotope;
    if(s->ranges[i_range+1].x - s->ranges[i_range].x == 0) /* Zero width. Return value of left side. */
        return s->cbins[i];
    return s->cbins[i] + ((s->cbins[i+s->n_isotopes] - s->cbins[i])/(s->ranges[i_range+1].x - s->ranges[i_range].x)) * (x - s->ranges[i_range].x);
}

sample *sample_alloc(size_t n_isotopes, size_t n_ranges) {
    sample *s = malloc(sizeof(sample));
    if(!s)
        return NULL;
    s->no_conc_gradients = FALSE;
    s->n_ranges = n_ranges;
    s->n_isotopes = n_isotopes;
    s->isotopes = calloc(n_isotopes, sizeof(jibal_isotope *));
    s->cbins = calloc( n_ranges * n_isotopes, sizeof(double));
    s->ranges = calloc(n_ranges, sizeof(struct sample_range));
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

int sample_model_sanity_check(const sample_model *sm) {
    for(size_t i = 0; i < sm->n_ranges; i++) {
        const sample_range *r = &sm->ranges[i];
        if(sm->type == SAMPLE_MODEL_LAYERED) {
            if(r->x < 0.0) {
                jabs_message(MSG_ERROR, stderr, "Sample model fails sanity check (layer %zu thickness negative).\n", i + 1);
                return EXIT_FAILURE;
            }
        } else if(sm->type == SAMPLE_MODEL_POINT_BY_POINT) {
            if(i && r->x < sm->ranges[i-1].x) {
                jabs_message(MSG_ERROR, stderr, "Sample model fails sanity check (non-monotonous range at %zu).\n", i + 1);
                return EXIT_FAILURE;
            }
        }
        for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
            if(*(sample_model_conc_bin(sm, i, i_mat)) < 0.0) {
                jabs_message(MSG_ERROR, stderr, "Sample model fails sanity check (negative concentration of %s (%zu) in layer/range %zu).   \n", sm->materials[i_mat]->name, i_mat + 1, i + 1);
                return EXIT_FAILURE;
            }
        }
        if(r->rough.x < 0.0) {
            jabs_message(MSG_ERROR, stderr, "Sample model fails sanity check (negative roughness range at %zu).\n", i + 1);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

void sample_model_renormalize(sample_model *sm) {
#ifdef SAMPLE_MODEL_RENORM_MATERIALS
    /* This step should be unnecessary, unless materials are modified after they are created. */
    for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
        jibal_material_normalize(sm->materials[i_mat]);
    }
#endif
    for(size_t i = 0; i < sm->n_ranges; i++) {
        //fprintf(stderr, "Range %zu: %g tfu\n", i, sm->ranges[i].x/C_TFU);
        if(sm->type == SAMPLE_MODEL_LAYERED) {
            if(sm->ranges[i].x < 0.0)
                sm->ranges[i].x = 0.0;
        } else if (sm->type == SAMPLE_MODEL_POINT_BY_POINT) {
            if(i && sm->ranges[i].x < sm->ranges[i-1].x)
                sm->ranges[i].x = sm->ranges[i-1].x; /* TODO: is this enough? */
        }
        double sum = 0.0;
        for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
            if(*(sample_model_conc_bin(sm, i, i_mat)) < 0.0)
                *(sample_model_conc_bin(sm, i, i_mat)) = 0.0;
            sum += *(sample_model_conc_bin(sm, i, i_mat));
        }
        if(sum == 0.0)
            continue;
        for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
            *(sample_model_conc_bin(sm, i, i_mat)) /= sum;
        }
    }
}

void sample_renormalize(sample *sample) {
    for(size_t i = 0; i < sample->n_ranges; i++) {
        double sum = 0.0;
        for(size_t i_isotope = 0; i_isotope < sample->n_isotopes; i_isotope++) {
            if(*(sample_conc_bin(sample, i, i_isotope)) < 0.0)
                *(sample_conc_bin(sample, i, i_isotope)) = 0.0;
            sum += *(sample_conc_bin(sample, i, i_isotope));
        }
        if(sum == 0.0)
            continue;
        for(size_t i_isotope = 0; i_isotope < sample->n_isotopes; i_isotope++) {
            *(sample_conc_bin(sample, i, i_isotope)) /= sum;
        }
    }
}

sample_model *sample_model_split_elements(const sample_model *sm) {
    if(!sm)
        return NULL;
    sample_model *out = sample_model_alloc(sample_model_element_count(sm), sm->n_ranges);
    out->type = sm->type;
#ifdef DEBUG
    fprintf(stderr, "Splitting sample model with %zu ranges and %zu materials to one with %zu materials (=elements).\n", sm->n_ranges, sm->n_materials, out->n_materials);
#endif
    size_t i = 0; /* Material index in output */
    for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
        for(size_t i_elem = 0; i_elem < sm->materials[i_mat]->n_elements; i_elem++) {
            jibal_element *e = jibal_element_copy(&sm->materials[i_mat]->elements[i_elem], JIBAL_ALL_ISOTOPES);
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
            sm_out->ranges[i].rough.model = ROUGHNESS_NONE;
            sm_out->ranges[i].rough.n = 0;
            if(i) {
                sm_out->ranges[i].x = sm_out->ranges[i-1].x; /* In layer models, the first point has the same depth as the second point of previous layer */
                sm_out->ranges[i+1].x += sm_out->ranges[i].x; /* And the second point is naturally deeper by thickness, which was stored as "depth" temporarily */
            } else {
                sm_out->ranges[i].x = 0.0;
            }
        }
        sm_out->type = SAMPLE_MODEL_LAYERED; /* Yes, this is odd. */
    }
    return sm_out;
}

sample *sample_from_sample_model(const sample_model *sm) { /* TODO: renormalize concentrations! */
    if(!sm)
        return NULL;
#ifdef DEBUG
    fprintf(stderr, "Sample model is type %i, it has %zu materials and %zu ranges.\n", sm->type, sm->n_materials, sm->n_ranges);
    for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
        fprintf(stderr, "Material %zu is %s. There are %zu elements.\n", i_mat, sm->materials[i_mat]->name, sm->materials[i_mat]->n_elements);
        for(size_t i_elem = 0; i_elem < sm->materials[i_mat]->n_elements; i_elem++) {
            fprintf(stderr, "  element %zu (%s: Z = %i): %zu isotopes.\n",  i_mat, sm->materials[i_mat]->elements[i_elem].name, sm->materials[i_mat]->elements[i_elem].Z, sm->materials[i_mat]->elements[i_elem].n_isotopes);
            for(size_t i_isotope = 0; i_isotope < sm->materials[i_mat]->elements[i_elem].n_isotopes; i_isotope++) {
                const jibal_isotope *isotope = sm->materials[i_mat]->elements[i_elem].isotopes[i_isotope];
                fprintf(stderr, "    isotope %zu: %s, A=%i, abundance %g, conc %g\n", i_isotope, isotope->name, isotope->A, isotope->abundance, sm->materials[i_mat]->elements[i_elem].concs[i_isotope]);
            }
        }
    }
#endif
    sample *s = malloc(sizeof(sample));
    sample_model *sm_copy = NULL;
    if(sm->type == SAMPLE_MODEL_LAYERED) {
        sm_copy = sample_model_to_point_by_point(sm);
        if(!sm_copy) {
            free(s);
            return NULL;
        }
        sm = sm_copy;
        s->no_conc_gradients = TRUE;
    } else {
        s->no_conc_gradients = FALSE;
    }
#ifdef DEBUG
    fprintf(stderr, "Point-by-point sample model for simulation:\n");
    sample_model_print(NULL, sm);
#endif
    size_t n_isotopes = 0;
    for(size_t i_mat = 0; i_mat < sm->n_materials; i_mat++) {
        for(size_t i_elem = 0; i_elem < sm->materials[i_mat]->n_elements; i_elem++) {
            n_isotopes += sm->materials[i_mat]->elements[i_elem].n_isotopes;
#ifdef DEBUG
            fprintf(stderr, "Material %zu element %zu has %zu isotopes.\n", i_mat, i_elem, sm->materials[i_mat]->elements[i_elem].n_isotopes);
#endif
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
    s->n_isotopes = n_isotopes;

    s->n_ranges = sm->n_ranges;
    s->isotopes = isotopes;
#ifdef DEBUG
    fprintf(stderr, "%zu isotopes before duplicates have been removed, they are:\n", s->n_isotopes);
    for(i = 0; i < s->n_isotopes; i++) {
        fprintf(stderr, "%zu: %s  (Z = %i, A = %i)\n", i, s->isotopes[i]->name, s->isotopes[i]->Z, s->isotopes[i]->A);
    }
#endif
    sample_sort_and_remove_duplicate_isotopes(s);
#ifdef DEBUG
    fprintf(stderr, "%zu isotopes remain after duplicates have been removed, they are:\n", s->n_isotopes);
    for(i = 0; i < s->n_isotopes; i++) {
        fprintf(stderr, "%zu: %s  (Z = %i, A = %i)\n", i, s->isotopes[i]->name, s->isotopes[i]->Z, s->isotopes[i]->A);
    }
#endif
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

    for(size_t i_range = 0; i_range < s->n_ranges; i_range++) { /* Set defaults for roughness */
        sample_range *r = &s->ranges[i_range];
        if(r->rough.x < ROUGH_TOLERANCE) {
            r->rough.model = ROUGHNESS_NONE;
            r->rough.x = 0.0;
        }
        if(r->rough.model == ROUGHNESS_GAMMA && r->rough.n == 0) { /* Zero is not valid number, but it means "auto" */
            /* TODO: implement variable number of spectra based on absolute and relative roughness. */
            r->rough.n = GAMMA_ROUGHNESS_STEPS;
#ifdef DEBUG
            fprintf(stderr, "Range %zu roughness is gamma, number of steps is automatic and set to %zu\n", i_range, r->rough.n);
#endif
        }
    }

#ifdef DEBUG
    sample_print(NULL, s, 0);
#endif
    free(sm_copy);
    return s;
}

int sample_model_print(const char *filename, const sample_model *sm) {
    if(!sm)
        return EXIT_FAILURE;
    FILE *f = fopen_file_or_stream(filename, "w");
    if(!f) {
        return EXIT_FAILURE;
    }
    size_t n_rl = sample_model_number_of_rough_ranges(sm);
    size_t n_layers_corrections = sample_model_number_of_ranges_with_non_unity_corrections(sm);
    size_t n_nonzero_density = sample_model_number_of_range_with_non_zero_density(sm);
    switch(sm->type) {
        case SAMPLE_MODEL_NONE:
            jabs_message(MSG_INFO, f, "Sample model is none.\n"); /* Not an error as such. No output (to f) is created. */
            return 0;
            break;
        case SAMPLE_MODEL_POINT_BY_POINT:
            jabs_message(MSG_INFO, f, "       depth");
            break;
        case SAMPLE_MODEL_LAYERED:
            jabs_message(MSG_INFO, f, "       thick");
            break;
    }
    if(n_rl) {
        jabs_message(MSG_INFO, f, "        rough");
        jabs_message(MSG_INFO, f, " n_rough");
    }
    if(n_nonzero_density) {
        jabs_message(MSG_INFO, f, "  density");
    }
    if(n_layers_corrections) {
        jabs_message(MSG_INFO, f, "  yield");
        jabs_message(MSG_INFO, f, "  bragg");
        jabs_message(MSG_INFO, f, "  stragg");
    }
    for (size_t i = 0; i < sm->n_materials; i++) {
        jabs_message(MSG_INFO, f, " %8s", sm->materials[i]->name);
    }
    jabs_message(MSG_INFO, f, "\n");
    for (size_t i = 0; i < sm->n_ranges; i++) {
        const sample_range *r = &(sm->ranges[i]);
        jabs_message(MSG_INFO, f, "%12.3lf", r->x/C_TFU);
        if(n_rl) {
            jabs_message(MSG_INFO, f, " %12.3lf", r->rough.x / C_TFU);
            jabs_message(MSG_INFO, f, " %7zu", r->rough.n);
        }
        if(n_nonzero_density) {
            jabs_message(MSG_INFO, f, " %8.4lf", r->density / C_G_CM3);
        }
        if(n_layers_corrections) {
            jabs_message(MSG_INFO, f, " %5.4lf %5.4lf %6.4lf", r->yield, r->bragg, r->stragg);
        }
        for (size_t j = 0; j < sm->n_materials; j++) {
            jabs_message(MSG_INFO, f, " %8.4lf", *sample_model_conc_bin(sm, i, j) * 100.0);
        }
        jabs_message(MSG_INFO, f, "\n");
    }
    fclose_file_or_stream(f);
    return 0;
}

size_t sample_model_number_of_rough_ranges(const sample_model *sm) {
    if(!sm)
        return 0;
    size_t n = 0;
    for(size_t i = 0; i < sm->n_ranges; i++) {
        if(sm->ranges[i].rough.model != ROUGHNESS_NONE)
            n++;
    }
    return n;
}

size_t sample_model_number_of_ranges_with_non_unity_corrections(const sample_model *sm) { /* Used to determine if "yield", "bragg" or "stragg" fields need to be output by sample_model_print() */
    if(!sm)
        return 0;
    size_t n = 0;
    for(size_t i = 0; i < sm->n_ranges; i++) {
        const sample_range *r = &(sm->ranges[i]);
        if(r->yield != 1.0 || r->bragg != 1.0 || r->stragg != 1.0)
            n++;
    }
    return n;
}

size_t sample_model_number_of_range_with_non_zero_density(const sample_model *sm) { /* Used to determine if "nm" field needs to be output by sample_model_print() */
    if(!sm)
        return 0;
    size_t n = 0;
    for(size_t i = 0; i < sm->n_ranges; i++) {
        const sample_range *r = &(sm->ranges[i]);
        if(r->density > 0.0)
            n++;
    }
    return n;
}

sample_model *sample_model_from_file(const jibal *jibal, const char *filename) {
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

    size_t i_depth = 0, i_rough = 0, i_nrough = 0, i_bragg = 0, i_yield = 0, i_stragg = 0, i_density = 0;

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
        while((col = strsep_with_quotes(&line_split, " \t")) != NULL) {
            if(*col == '\0') {/* Multiple separators are treated as one */
                continue;
            }
            if(headers) { /* Headers */
                if(strncmp(col, "bragg", 5) == 0) {
                    i_bragg = n;
                } else if(strncmp(col, "yield", 5) == 0) {
                    i_yield = n;
                } else if(strncmp(col, "stragg", 5) == 0) {
                    i_stragg = n;
                } else if(strncmp(col, "density", 7) == 0) {
                    i_density = n;
                } else if(strncmp(col, "rough", 5) == 0) {
                    i_rough = n;
                } else if(strncmp(col, "n_rough", 7) == 0) {
                    i_nrough = n;
                } else if(strncmp(col, "thick", 5) == 0) {
                    i_depth = n;
                    sm->type = SAMPLE_MODEL_LAYERED;
                } else if(strncmp(col, "depth", 5) == 0) {
                    i_depth = n;
                    sm->type = SAMPLE_MODEL_POINT_BY_POINT;
                } else {
                    sm->materials = realloc(sm->materials, sizeof (jibal_material *) * (sm->n_materials + 1)); /* TODO: check allocation */
                    sm->materials[sm->n_materials] = jibal_material_create(jibal->elements, col);
                    if(!sm->materials[sm->n_materials]) {
                        jabs_message(MSG_WARNING, stderr, "Could not create material %s, ignoring it!\n", col);
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
                sample_range *r = &sm->ranges[sm->n_ranges];
                r->x = 0.0;
                r->yield = 1.0;
                r->bragg = 1.0;
                r->stragg = 1.0;
                r->density = 0.0;
                r->rough.x = 0.0;
                r->rough.model = ROUGHNESS_NONE;
                r->rough.n = 0;
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
            } else if (n == i_bragg) {
                r->bragg = x;
            } else if (n == i_yield) {
                r->yield = x;
            } else if (n == i_stragg) {
                r->stragg = x;
            } else if (n == i_rough) {
                r->density = x*C_G_CM3;
            } else if (n == i_rough) {
                r->rough.x = x*C_TFU;
            } else if (n == i_nrough) {
                r->rough.n = floor(x);
            }else if(i_material < sm->n_materials) {
                *(sample_model_conc_bin(sm, sm->n_ranges-1, i_material)) = x;
                i_material++;
            }
            n++;
        }
        headers = 0;
    }

    sample_model_renormalize(sm);

    for(size_t i_range = 0; i_range < sm->n_ranges; i_range++) { /* Set defaults for roughness*/
        sample_range *r = &sm->ranges[i_range];
        if(r->rough.x > 0.1*C_TFU) {
            r->rough.model = ROUGHNESS_GAMMA; /* TODO: other models */
        } else {
            r->rough.model = ROUGHNESS_NONE;
            r->rough.x = 0.0;
        }
    }
    free(line);
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

sample_model *sample_model_from_argv(const jibal *jibal, int * const argc, char * const ** const argv) {
    sample_model *sm = malloc(sizeof(sample_model));
    sm->type = SAMPLE_MODEL_LAYERED;
    size_t n = 0;
    sm->n_ranges = 0;
    sm->n_materials = 0;
    sm->materials = NULL;
    sm->ranges = NULL;
    sm->cbins = NULL;
    int simplify = TRUE;
    if((*argc) >= 1) {
        if(strcmp((*argv)[0], "nosimplify") == 0) {
            simplify = FALSE;
            (*argc)--;
            (*argv)++;
        }
    }
    while ((*argc) >= 2) {
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
        if(sm->n_ranges && strcmp((*argv)[0], "rough") == 0) {
            sample_range *range = &sm->ranges[sm->n_ranges - 1];
            range->rough.x = jibal_get_val(jibal->units, UNIT_TYPE_LAYER_THICKNESS, (*argv)[1]);
            range->rough.model = ROUGHNESS_GAMMA;
        } else if(sm->n_ranges && strcmp((*argv)[0], "n_rough") == 0) {
                sample_range *range = &sm->ranges[sm->n_ranges - 1];
                range->rough.n = strtoul((*argv)[1], NULL, 10);
        } else if(sm->n_ranges && strcmp((*argv)[0], "yield") == 0) {
                sample_range *range = &sm->ranges[sm->n_ranges - 1];
                range->yield = strtod((*argv)[1], NULL);
        } else if(sm->n_ranges && strcmp((*argv)[0], "bragg") == 0) {
                sample_range *range = &sm->ranges[sm->n_ranges - 1];
                range->bragg = strtod((*argv)[1], NULL);
        } else if(sm->n_ranges && strcmp((*argv)[0], "stragg") == 0) {
            sample_range *range = &sm->ranges[sm->n_ranges - 1];
            range->stragg = strtod((*argv)[1], NULL);
        } else if(sm->n_ranges && strcmp((*argv)[0], "density") == 0) {
            sample_range *range = &sm->ranges[sm->n_ranges - 1];
            range->density = jibal_get_val(jibal->units, UNIT_TYPE_DENSITY, (*argv)[1]);
        } else {
            sm->materials[sm->n_ranges] = jibal_material_create(jibal->elements, (*argv)[0]);
            if(!sm->materials[sm->n_ranges]) {
#ifdef DEBUG
                fprintf(stderr, "Material from formula \"%s\" was NOT created. Finishing after %zu ranges and %zu materials.\n", (*argv)[0], sm->n_ranges, sm->n_materials);
#endif
                break;
            }
#ifdef DEBUG
            fprintf(stderr, "Material from formula \"%s\" was created\n", (*argv)[0]);
#endif
            sample_range *range = &sm->ranges[sm->n_ranges];
            range->x = jibal_get_val(jibal->units, UNIT_TYPE_LAYER_THICKNESS, (*argv)[1]);
            range->bragg = 1.0;
            range->yield = 1.0;
            range->stragg = 1.0;
            range->density = 0.0;
            range->rough.x = 0.0;
            range->rough.model = ROUGHNESS_NONE;
            range->rough.n = 0;
            sm->n_ranges++;
            sm->n_materials++;
        }
        (*argc) -= 2;
        (*argv) += 2;
    }

    if(sm->n_ranges == 0) {
#ifdef DEBUG
        fprintf(stderr, "No ranges were parsed.\n");
#endif
        sample_model_free(sm);
        return NULL;
    }

    sm->cbins = calloc(sm->n_ranges * sm->n_ranges, sizeof(double));
    for(size_t i = 0; i < sm->n_ranges; i++) {
        *sample_model_conc_bin(sm, i, i) = 1.0;
    }
    if(simplify) {
#ifdef DEBUG
        fprintf(stderr, "Sample model from argv before splitting elements:\n");
        sample_model_print(NULL, sm);
#endif
        sample_model *sm2 = sample_model_split_elements(sm);
        sample_model_free(sm);
        sm = sm2;
#ifdef DEBUG
        fprintf(stderr, "Sample model after splitting elements:\n");
        sample_model_print(NULL, sm2);
#endif
    }
    return sm;
}

sample_model *sample_model_from_string(const jibal *jibal, const char *str) {
    int argc_orig = 0;
    char **argv_orig = string_to_argv(str, &argc_orig);
    char **argv = argv_orig;
    int argc = argc_orig;
    sample_model *sm = sample_model_from_argv(jibal, &argc, (char *const **) &argv);
#ifdef DEBUG
    fprintf(stderr, "%i arguments remain after sample conversion\n", argc);
#endif
    argv_free(argv_orig, argc_orig);
    return sm;
}

char *sample_model_to_string(const sample_model *sm) {
    char *out = NULL;
    if(sm->type != SAMPLE_MODEL_LAYERED) /* TODO: point-by-point */
        return NULL;
    for(size_t i_range= 0; i_range < sm->n_ranges; i_range++) {
        const sample_range *r = &sm->ranges[i_range];
        asprintf_append(&out, "%s", i_range?" ":"");
        for(size_t i_material = 0; i_material < sm->n_materials; i_material++) {
            double conc = *(sample_model_conc_bin(sm, i_range, i_material));
            const char *name = sm->materials[i_material]->name;
            if(conc > (1.0 - CONC_TOLERANCE)) { /* Omit unnecessary "1" */
                asprintf_append(&out, "%s", name);
            } else if(conc > CONC_TOLERANCE) { /* Only print if concentration above zero. */
                asprintf_append(&out, "%s%g", name, conc);
            }
        }
        asprintf_append(&out, " %g%s", r->x/C_TFU, "tfu");
        if(r->rough.model != ROUGHNESS_NONE && r->x > ROUGH_TOLERANCE) {
            asprintf_append(&out, " rough %gtfu", r->rough.x/C_TFU, "tfu");
            if(r->rough.n != 0) {
                asprintf_append(&out, " n_rough %zu", r->rough.n);
            }
        }
        if(r->bragg != 1.0) {
            asprintf_append(&out, " bragg %g", r->bragg);
        }
        if(r->yield != 1.0) {
            asprintf_append(&out, " yield %g", r->yield);
        }
        if(r->stragg != 1.0) {
            asprintf_append(&out, " stragg %g", r->stragg);
        }
        if(r->density != 0.0) {
            asprintf_append(&out, " density %g", r->density);
        }
    }
    return out;
}

sample *sample_copy(const sample *s_in) {
    sample *s_out = sample_alloc(s_in->n_isotopes, s_in->n_ranges);
    if(!s_out)
        return NULL;
    s_out->no_conc_gradients = s_in->no_conc_gradients;
    memcpy(s_out->isotopes, s_in->isotopes, sizeof(jibal_isotope *) * s_out->n_isotopes);
    memcpy(s_out->ranges, s_in->ranges, sizeof (struct sample_range) * s_out->n_ranges);
    memcpy(s_out->cbins, s_in->cbins, sizeof (double) * s_out->n_isotopes * s_out->n_ranges);
#ifdef DEBUG
        fprintf(stderr, "Made a sample copy with %zu ranges and %zu isotopes.\n", s_out->n_ranges, s_out->n_isotopes);
#endif
    return s_out;
}

double sample_mass_density_range(const sample *sample, size_t i_range) {
    if(i_range < 1 || i_range >= sample->n_ranges)
        return 0.0;
    double thickness = (sample->ranges[i_range].x - sample->ranges[i_range-1].x); /* areal density */
    double sum = 0.0;
    for (size_t i = 0; i < sample->n_isotopes; i++) {
        double conc = 0.5 * ( *(sample_conc_bin(sample, i_range, i)) + *(sample_conc_bin(sample, i_range-1, i))); /* average */
        sum += conc * sample->isotopes[i]->mass;
    }
    return sum * thickness;
}

double sample_thickness_in_nm_range(const sample *sample, size_t i_range) {
    if(i_range < 1 || i_range >= sample->n_ranges)
        return 0.0;
    double massdensity = sample_mass_density_range(sample, i_range); /* kg/m2 areal mass density */
    double avgdensity = 0.5*(sample->ranges[i_range-1].density + sample->ranges[i_range].density); /* kg/m3 density */
    return massdensity/avgdensity;
}

double sample_areal_density_isotope_range(const sample *sample, size_t i_isotope, size_t i_range) {
    if(i_range < 1 || i_range >= sample->n_ranges)
        return 0.0;
    if(i_isotope >= sample->n_isotopes)
        return 0.0;
    double thickness = (sample->ranges[i_range].x - sample->ranges[i_range-1].x);
    return 0.5 * ( *(sample_conc_bin(sample, i_range, i_isotope)) + *(sample_conc_bin(sample, i_range-1, i_isotope))) * thickness;
}

void sample_areal_densities_print(FILE *f, const sample *sample, int print_isotopes) {
    if(!sample)
        return;
    jabs_message(MSG_INFO, f, "SUM (tfu)        ");
    double sum = 0.0;
    for (size_t i = 0; i < sample->n_isotopes; i++) {
        for (size_t j = 1; j < sample->n_ranges; j++) {
            sum += sample_areal_density_isotope_range(sample, i, j);
        }
        if (print_isotopes || i == sample->n_isotopes-1 || sample->isotopes[i]->Z != sample->isotopes[i+1]->Z) {
            jabs_message(MSG_INFO, f, " %9.3lf", sum/C_TFU);
            sum = 0.0;
        }
    }
    jabs_message(MSG_INFO, f, "\n");
}


int sample_print_thicknesses(const char *filename, const sample *sample) {
    if(!sample)
        return EXIT_FAILURE;
    FILE *f = fopen_file_or_stream(filename, "w");
    if(!f)
        return EXIT_FAILURE;
    if(sample->no_conc_gradients) {
        jabs_message(MSG_INFO, f, "Layer       tfu     ug/cm2       nm\n");
        for(size_t i = 1; i < sample->n_ranges; i += 2) {
            jabs_message(MSG_INFO, f, "%5zu %9.3lf  %9.3lf %8.1lf\n", (i + 1)/2,
                         (sample->ranges[i].x - sample->ranges[i-1].x) / C_TFU,
                         sample_mass_density_range(sample, i) / (C_UG / C_CM2),
                         sample_thickness_in_nm_range(sample, i) / C_NM);
        }
    }
    fclose_file_or_stream(f);
    return EXIT_SUCCESS;
}

int sample_print(const char *filename, const sample *sample, int print_isotopes) {
    if(!sample)
        return EXIT_FAILURE;
    FILE *f = fopen_file_or_stream(filename, "w");
    if(!f)
        return EXIT_FAILURE;
    jabs_message(MSG_INFO, f, "    DEPTH   ROUGH");
    int Z = 0;
    for (size_t i = 0; i < sample->n_isotopes; i++) {
        if(print_isotopes) {
            jabs_message(MSG_INFO, f, " %9s", sample->isotopes[i]->name);
        } else if(Z != sample->isotopes[i]->Z){
            const char *s = sample->isotopes[i]->name;
            while(*s >= '0' && *s <= '9') {s++;} /* Skip numbers, e.g. 28Si -> Si */
            jabs_message(MSG_INFO, f, " %9s", s);
            Z = sample->isotopes[i]->Z; /* New element */
        }
    }
    jabs_message(MSG_INFO, f, "\n");
    jabs_message(MSG_INFO, f, "      tfu     tfu\n");
    for (size_t i = 0; i < sample->n_ranges; i++) {
        jabs_message(MSG_INFO, f, "%9.3lf", sample->ranges[i].x/C_TFU);
        jabs_message(MSG_INFO, f, " %7.3lf", sample->ranges[i].rough.x/C_TFU);
        double sum = 0.0;
        for (size_t j = 0; j < sample->n_isotopes; j++) {
            sum += sample->cbins[i * sample->n_isotopes + j];
            if (print_isotopes || j == sample->n_isotopes-1 || sample->isotopes[j]->Z != sample->isotopes[j+1]->Z) { /* Last isotope or next isotope belongs to another element, print. */
                jabs_message(MSG_INFO, f, " %9.4lf", sum * 100.0);
                sum = 0.0;
            }
        }
        jabs_message(MSG_INFO, f, "\n");
    }
    jabs_message(MSG_INFO, f, "\n");
    fclose_file_or_stream(f);
    return EXIT_SUCCESS;
}

void sample_free(sample *sample) {
    if(!sample)
        return;
    free(sample->ranges);
    free(sample->isotopes);
    free(sample->cbins);
    free(sample);
}

double sample_isotope_max_depth(const sample *sample, size_t i_isotope) {
    for(size_t i = sample->n_ranges; i--;) {
        double *c = sample_conc_bin(sample, i, i_isotope);
        if(*c > ABUNDANCE_THRESHOLD) {
            return sample->ranges[i].x;
        }
    }
    return 0.0;
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
