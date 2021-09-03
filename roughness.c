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



thick_prob_dist *thickness_probability_table_gen(double thickness, double sigma, size_t n) {
    if(thickness < 0.0) {
        thickness = 0.0;
    }
    if(sigma < 0.01 * C_TFU) { /* Special case for (near) zero sigma, reduce n to 1 */
        sigma = 0.0;
        n = 1;
    }

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
    double areal_sum = 0.0;
    if(n == 1) {
        tpd->p[0].prob = 1.0;
        tpd->p[0].x = thickness;
        return tpd;
    }
    for(i = 0; i < n; i++) {
        thick_prob *p = &tpd->p[i];
        double x_low = low + i*step;
        double x_high = low + (i+1)*step;
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
    for(i = 0; i < n; i++) {
        thick_prob *p = &tpd->p[i];
        p->prob *= corr;
    }
    return tpd;
}

void thickness_probability_table_free(thick_prob_dist *tpd) {
    if(!tpd)
        return;
    free(tpd->p);
    free(tpd);
}
