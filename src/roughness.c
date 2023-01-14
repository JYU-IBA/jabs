/*

    Jaakko's Backscattering Simulator (JaBS)
    Copyright (C) 2021 - 2023 Jaakko Julin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    See LICENSE.txt for the full license.

 */

#include <stdio.h>
#include <string.h>
#include <jibal_units.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "generic.h"
#include "defaults.h"
#include "roughness.h"

double thickness_gamma_pdf(double x, double thickness, double sigma) {
    const double a = (thickness*thickness)/(sigma*sigma);
    const double b = (sigma*sigma)/thickness;

    return gsl_cdf_gamma_P(x, a, b);
}

thick_prob_dist *thickness_probability_table_gamma(double thickness, double sigma, size_t n) {
    if(thickness < 0.0) {
        thickness = 0.0;
    }
    if(sigma < 0.01 * C_TFU) { /* Special case for (near) zero sigma, reduce n to 1 */
        sigma = 0.0;
        n = 1;
    }

    double low = thickness - sigma*4.0;
    double high = thickness + sigma*4.0;
    if(low < 0.0) {
        low = 0.0;
    }
    double step = (high-low)/(n*1.0);
    thick_prob_dist *tpd = thickness_probability_table_new(n);
    double sum = 0.0;
    double areal_sum = 0.0;
    if(n == 1) {
        tpd->p[0].prob = 1.0;
        tpd->p[0].x = thickness;
        return tpd;
    }
    for(size_t i = 0; i < n; i++) {
        thick_prob *p = &tpd->p[i];
        double x_low = low + i*step;
        double x_high = x_low + step;
        p->x = (x_high + x_low) / 2.0;
        p->prob = thickness_gamma_pdf(x_high, thickness, sigma);
        p->prob -= thickness_gamma_pdf(x_low, thickness, sigma);
        sum += p->prob;
        areal_sum += p->prob * p->x;
    }
    double corr = thickness/areal_sum; /* Correction factor, since the used distribution is not perfect due to cutoffs. */
#ifdef DEBUG
    fprintf(stderr, "Gamma roughness, thickness low %g, high %g, step %g. Sum of probabilities %g, areal density weighted with probability is %g tfu. Correction factor %g will be applied.\n", low/C_TFU, high/C_TFU, step/C_TFU, sum, areal_sum/C_TFU, corr);
#endif
    for(size_t i = 0; i < n; i++) {
        thick_prob *p = &tpd->p[i];
        p->prob *= corr;
    }
    return tpd;
}


thick_prob_dist *thickness_probability_table_new(size_t n) {
    thick_prob_dist *tpd = malloc(sizeof(thick_prob_dist));
    tpd->modulo = 0;
    tpd->i_range = 0;
    tpd->n = n;
    tpd->p = calloc(n, sizeof(thick_prob));
    return tpd;
}

thick_prob_dist *thickness_probability_table_from_file(const char *filename) { /* TODO: untested, contains bugs with extremely high probability! */
    size_t n = ROUGHNESS_SUBSPECTRA_MAXIMUM, n_true = 0;
    char *line = NULL;
    size_t line_size = 0;
    size_t lineno = 0;
    int fail = FALSE;
    FILE *in;
    if(!filename)
        return NULL;
    in = fopen(filename, "r");
    if(!in)
        return NULL;
    thick_prob_dist *tpd = thickness_probability_table_new(n);
    while(getline(&line, &line_size, in) > 0) {
        lineno++;
        line[strcspn(line, "\r\n")] = 0; /* Strips all kinds of newlines! */
        if(strlen(line) >= 1 && *line == '#') {/* Comment */
            continue;
        }
        if(n_true == n) {
            fail = TRUE;
#ifdef DEBUG
            fprintf(stderr, "TPD file too long, max = %zu.\n", n);
#endif
            break;
        }
        char *s = line, *end;
        double thick = strtod(s, &end);
        if(s == end) {
            fail = TRUE;
            break;
        }
        s = end;
        double prob = strtod(s, &end);
        if(*end != '\0') {
            fail = TRUE;
            break;
        }
        thick_prob *p = &(tpd->p[n_true]);
        p->prob = prob;
        p->x = thick * C_TFU;
#ifdef DEBUG
        fprintf(stderr, "Line %zu: parsed \"%s\" into %g tfu and %g\n", lineno, line, p->x/C_TFU, p->prob);
#endif
        n_true++;
    }
    fclose(in);
    if(n_true == 0) {
        fail = TRUE;
    }
    if(fail) {
#ifdef DEBUG
        fprintf(stderr, "Failure in reading line %zu of file \"%s\".\n", lineno, filename);
#endif
        thickness_probability_table_free(tpd);
        return NULL;
    }
    thickness_probability_table_realloc(tpd, n_true); /* Shrink to true size */
    return tpd;
}

void thickness_probability_table_normalize(thick_prob_dist *tpd) {
    if(!tpd)
        return;
    double psum = 0.0;
    for(size_t i = 0; i < tpd->n; i++) {
        thick_prob *p = &tpd->p[i];
        psum += p->prob;
    }
#ifdef DEBUG
    fprintf(stderr, "TPD sum of probabilities is %g.\n", psum);
#endif
    for(size_t i = 0; i < tpd->n; i++) {
        thick_prob *p = &tpd->p[i];
        p->prob /= psum;
    }
}

double thickness_probability_table_areal_density(thick_prob_dist *tpd) { /* Doesn't make sense unless the table is normalized */
    if(!tpd)
        return 0.0;
    double ad = 0.0;
    for(size_t i = 0; i < tpd->n; i++) {
        thick_prob *p = &tpd->p[i];
        ad += p->prob * p->x;
    }
    return ad;
}

void thickness_probability_table_free(thick_prob_dist *tpd) {
    if(!tpd)
        return;
    free(tpd->p);
    free(tpd);
}

void thickness_probability_table_print(FILE *f, const thick_prob_dist *tpd) {
    if(!tpd)
        return;
    fprintf(f, " i | thickness | weight(%%)\n");
    for(size_t i = 0; i < tpd->n; i++) {
        thick_prob *p = &tpd->p[i];
        fprintf(f, "%2zu | %9.2lf | %9.4lf\n", i, p->x/C_TFU, p->prob/C_PERCENT);
    }
}

thick_prob_dist *thickness_probability_table_realloc(thick_prob_dist *tpd, size_t n) {
    tpd->p = realloc(tpd->p, n * sizeof(thick_prob));
    if(!tpd->p) {
        thickness_probability_table_free(tpd);
        return NULL;
    }
    tpd->n = n;
    return tpd;
}

thick_prob_dist *thickness_probability_table_copy(const thick_prob_dist *tpd) {
    if(!tpd)
        return NULL;
    thick_prob_dist *tpd_copy = thickness_probability_table_new(tpd->n);
    for(size_t i = 0; i < tpd_copy -> n; i++) {
        tpd_copy->p[i] = tpd->p[i];
    }
    return tpd_copy;
}

int roughness_reset(roughness *r) {
    if(r->model == ROUGHNESS_FILE) {
        roughness_file_free(r->file);
    }
    r->model = ROUGHNESS_NONE;
    r->x = 0.0;
    r->n = 0;
    r->file = NULL;
    return EXIT_SUCCESS;
}

int roughness_reset_if_below_tolerance(roughness *r) {
    if(r->model == ROUGHNESS_FILE) { /* Files are kept as-is. */
        return EXIT_SUCCESS;
    }
    if(r->x < ROUGH_TOLERANCE) {
        return roughness_reset(r);
    }
    return EXIT_SUCCESS;
}

int roughness_set_from_file(roughness *r, const char *filename) {
    roughness_file *rf = roughness_file_read(filename);
    if(!rf) {
        return EXIT_FAILURE;
    }
    r->file = rf;
    r->model = ROUGHNESS_FILE;
    r->n = rf->tpd->n; /* "n" should always match, but in calculations tpd->n takes precedence. r->n is for convenience, since it is no longer a variable for a roughness model. */
    r->x = 0.0;
    return EXIT_SUCCESS;
}

roughness_file *roughness_file_read(const char *filename) {
    if(!filename)
        return NULL;
    thick_prob_dist *tpd = thickness_probability_table_from_file(filename);
    if(!tpd)
        return NULL;
    roughness_file *rf = malloc(sizeof(roughness_file));
    rf->filename = strdup(filename);
    rf->tpd = tpd;
    thickness_probability_table_normalize(rf->tpd);
#ifdef DEBUG
    thickness_probability_table_print(stderr, rf->tpd);
#endif
    return rf;
}

roughness_file *roughness_file_copy(const roughness_file *rf) {
    if(!rf)
        return NULL;
    roughness_file *rf_out = malloc(sizeof(roughness_file));
    rf_out->tpd = thickness_probability_table_copy(rf->tpd);
    rf_out->filename = strdup_non_null(rf->filename);
    return rf_out;
}

void roughness_file_free(roughness_file *rf) {
    if(!rf)
        return;
    free(rf->filename);
    thickness_probability_table_free(rf->tpd);
    free(rf);
}
