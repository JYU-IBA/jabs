#include <math.h>
#include <assert.h>
#include <jibal_units.h>
#include <gsl/gsl_sf_erf.h>
#include "brick.h"

extern inline double erf_Q_fast(double x) { /* Approximative gaussian CDF */
    return x < 0.0 ? 1.0-0.5*exp(0.77428768622*x-0.37825569191*x*x) : 0.5*exp(-0.77428768622*x-0.37825569191*x*x);
    /* Based on: Tsay, WJ., Huang, C.J., Fu, TT. et al. J Prod Anal 39, 259–269 (2013). https://doi.org/10.1007/s11123-012-0283-1 */
}

void bricks_calculate_sigma(const detector *det, const jibal_isotope *isotope, brick *bricks, size_t last_brick) {
    for(size_t i = 0; i <= last_brick; i++) {
        bricks[i].sigma = sqrt(bricks[i].S + detector_resolution(det, isotope, bricks[i].E) +  bricks[i].S_geo_x + bricks[i].S_geo_y);
    }
}

void bricks_convolute(gsl_histogram *h, const brick *bricks, size_t last_brick, const double scale, const double sigmas_cutoff, int accurate) {
    double (*erf_Q)(double);
    if(accurate) {
        erf_Q = gsl_sf_erf_Q;
    } else {
        erf_Q = erf_Q_fast;
    }

    for(size_t i = 1; i <= last_brick; i++) {
        const brick *b_high = &bricks[i-1];
        const brick *b_low = &bricks[i];
        double E_cutoff_low = b_low->E - b_low->sigma * sigmas_cutoff; /* Low energy cutoff (brick) */
        double E_cutoff_high = b_high->E + b_high->sigma * sigmas_cutoff;
        double b_w_inv = 1.0/(b_high->E - b_low->E); /* inverse of brick width (in energy) */
        size_t lo = 0, mi, hi = h->n;
        while(hi - lo > 1) { /* Find histogram range with E_cutoff. This should be safe. */
            mi = (hi + lo) / 2;
            if(E_cutoff_low >= h->range[mi]) {
                lo = mi;
            } else {
                hi = mi;
            }
        }
        for(size_t j = lo; j < h->n; j++) {
            if(h->range[j] > E_cutoff_high) /* Low energy edge of histogram is above cutoff */
                break; /* Assumes histograms have increasing energy */
            const double E = (h->range[j] + h->range[j + 1]) / 2.0; /* Approximate gaussian at center bin */
            const double w = h->range[j + 1] - h->range[j];
            const double y = (erf_Q((b_low->E - E) / b_low->sigma) - erf_Q((b_high->E - E) / b_high->sigma));
            h->bin[j] += scale * y * w * b_w_inv * b_low->Q;
        }
    }
}
