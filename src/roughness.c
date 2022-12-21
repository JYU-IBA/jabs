#include <stdio.h>
#include <jibal_units.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <generic.h>
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

void thickness_probability_table_free(thick_prob_dist *tpd) {
    if(!tpd)
        return;
    free(tpd->p);
    free(tpd);
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
}
