#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "prob_dist.h"


prob_dist *prob_dist_alloc(size_t n) {
    prob_dist *pd = malloc(sizeof(prob_dist));
    if(!pd)
        return NULL;
    pd->n = n;
    pd->points = calloc(n, sizeof(prob_point));
    if(!pd->points) {
        free(pd);
        return NULL;
    }
    return pd;
}

void prob_dist_free(prob_dist *pd) {
    if(!pd)
        return;
    free(pd->points);
    free(pd);
}

prob_dist *prob_dist_gaussian(size_t n) { /* creates a discrete gaussian (mu = 0, sigma = 1) probability distribution, sum of probabilities is 100% */
    if(n < 1)
        return NULL;
    prob_dist *pd = prob_dist_alloc(n);
    if(!pd)
        return NULL;
    double sigmas = 0.5 * n; /* number of standard deviations (+/-) around mean. n dependent scaling here means the bin width is one standard deviation */
    if(sigmas >= 3.0) /* +/- 3.0 sigmas is large enough, use extra bins for extra accuracy */
        sigmas = 3.0;
    double low = -1.0 * sigmas;
    double high = 1.0 * sigmas;
    double step = (high - low)/n;
    double p_sum = 0.0;
    for(size_t i = 0; i < n; i++) {
        prob_point *pp = &(pd->points[i]);
        double x_low = low + i * step;
        double x_high = low + (i + 1) * step;
        double x = (x_low + x_high)/2.0;
        double p = gsl_cdf_gaussian_P(x_high, 1.0) - gsl_cdf_gaussian_P(x_low, 1.0);
        pp->x = x;
        pp->p = p;
#ifdef DEBUG_CS_VERBOSE
        fprintf(stderr, "%zu %12g %12g %12g %12g\n", i, pp->x, pp->p, x_low, x_high);
#endif
        p_sum += p;
    }
#ifdef DEBUG_CS_VERBOSE
    fprintf(stderr, "P_sum = %.7lf (will be compensated for)\n", p_sum);
#endif
    for(size_t i = 0; i < n; i++) {
        prob_point *pp = &(pd->points[i]);
        pp->p /= p_sum;
#ifdef DEBUG_CS_VERBOSE
        fprintf(stderr, "%zu %12g %12g\n", i, pp->x, pp->p);
#endif
    }
    return pd;
}
