#include <stdio.h>
#include <jibal_units.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "roughness.h"

double thickness_gamma_pdf(double x, double thickness, double sigma) {
    const double a = (thickness*thickness)/(sigma*sigma);
    const double b = (sigma*sigma)/thickness;

    return gsl_cdf_gamma_P(x, a, b);
    //return gsl_ran_gamma_pdf(x, a, b);
}



thick_prob_dist *thickness_probability_table_gen(double thickness, double sigma, const size_t n) {
    double low = thickness - sigma*4.0;
    double high = thickness + sigma*4.0;
    if(low < 0.0)
        low = 0.0;
    size_t i;
    double step = (high-low)/(n*1.0);
    thick_prob_dist *tpd = malloc(sizeof(thick_prob_dist));
    tpd->n = n;
    tpd->p = malloc(sizeof(thick_prob) * tpd->n);

    
    double sum = 0.0;
    for(i = 0; i < n; i++) {
        thick_prob *p = &tpd->p[i];
        double x_low = low + i*step;
        double x_high = low + (i+1)*step;
        p->x = (x_high + x_low) / 2.0;
        p->prob = thickness_gamma_pdf(x_high, thickness, sigma);
        p->prob -= thickness_gamma_pdf(x_low, thickness, sigma);
        sum += p->prob;
    }
#ifdef DEBUG
    fprintf(stderr, "Gamma roughness, thickness low %g, high %g, step %g. Sum of probabilities before normalization %g.\n", low/C_TFU, high/C_TFU, step/C_TFU, sum);
#endif
    for(i = 0; i < n; i++) {
        thick_prob *p = &tpd->p[i];
        p->prob /= sum; /* Normalize */
    }
    return tpd;
}

void thickness_probability_table_free(thick_prob_dist *tpd) {
    if(!tpd)
        return;
    free(tpd->p);
    free(tpd);
}
