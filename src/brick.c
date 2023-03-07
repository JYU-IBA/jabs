#include <math.h>
#include <assert.h>
#include <jibal_units.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_math.h>
#include "brick.h"

extern inline double erf_Q_fast(double x) { /* Approximative gaussian CDF */
    return x < 0.0 ? 1.0-0.5*exp(0.77428768622*x-0.37825569191*x*x) : 0.5*exp(-0.77428768622*x-0.37825569191*x*x);
    /* Based on: Tsay, WJ., Huang, C.J., Fu, TT. et al. J Prod Anal 39, 259–269 (2013). https://doi.org/10.1007/s11123-012-0283-1 */
}

void bricks_calculate_sigma(const detector *det, const jibal_isotope *isotope, brick *bricks, size_t last_brick) {
    for(size_t i = 0; i <= last_brick; i++) {
        //double old = bricks[i].S_sum;
        bricks[i].S_sum = sqrt(bricks[i].S + detector_resolution(det, isotope, bricks[i].E) + bricks[i].S_geo_x + bricks[i].S_geo_y);
        //fprintf(stderr, " %zu: old S_sum = %.12g keV, new = %.12g keV\n", i, old / C_KEV, bricks[i].S_sum / C_KEV);
    }
}

void bricks_convolute(jabs_histogram *h, const calibration *c, const brick *bricks, size_t last_brick, const double scale, const double sigmas_cutoff, double emin, int accurate) {
    double (*erf_Q)(double);
    if(accurate) {
        erf_Q = gsl_sf_erf_Q;
    } else {
        erf_Q = erf_Q_fast;
    }

    for(size_t i = 1; i <= last_brick; i++) {
        const brick *b_low = &bricks[i]; /* Low refers to lower index number (i), not that b_low->E < b_high->E! (incident ion has lower energy as function of i, but not necessarily the detected energy of the reaction product) */
        if(b_low->Q == 0.0) {
            continue;
        }
        const brick *b_high = &bricks[i-1];
        double E_low, E_cutoff_low, E_cutoff_high; /* Lower energy edge of brick, Low energy cutoff (gaussian), High energy cutoff (gaussian) */
        if(b_low->E < b_high->E) { /* Detected energy increasing as brick number increases */
            E_low = b_low->E;
            E_cutoff_low = b_low->E - b_low->S_sum * sigmas_cutoff;
            E_cutoff_high = b_high->E + b_high->S_sum * sigmas_cutoff;
        } else { /* Energy decreasing as brick number increases */
            E_low = b_high->E;
            E_cutoff_low = b_high->E - b_high->S_sum * sigmas_cutoff;
            E_cutoff_high = b_low->E + b_low->S_sum * sigmas_cutoff;
        }
        if(E_low < emin) { /* Brick is already partially below emin (without convolution), it's low enough! */
            continue;
        }
        E_cutoff_low = GSL_MAX_DBL(E_cutoff_low, emin);
        double b_w_inv = 1.0/(b_high->E - b_low->E); /* inverse of brick width (in energy) */
        size_t ch_start = calibration_inverse(c, E_cutoff_low, h->n);
#ifdef DEBUG_BRICK_OUTPUT
        double norm = 0.0;
#endif
        for(size_t j = ch_start; j < h->n - 1; j++) {
            if(h->range[j] > E_cutoff_high) /* Low energy edge of histogram is above cutoff */
                break; /* Assumes histograms have increasing energy */
            const double E = (h->range[j] + h->range[j + 1]) / 2.0; /* Approximate gaussian at center bin */
            const double w = h->range[j + 1] - h->range[j];
            const double y = (erf_Q((b_low->E - E) / b_low->S_sum) - erf_Q((b_high->E - E) / b_high->S_sum));
            double out = scale * y * w * b_w_inv * b_low->Q;
            h->bin[j] += out;
#ifdef DEBUG_BRICK_OUTPUT
            fprintf(stdout, "BRICK %3zu %4zu %9.3lf %9.3lf %12e %12e %12e %12e %12e %12e\n", i, j, h->range[j] / C_KEV, h->range[j + 1] / C_KEV, out, h->bin[j], scale, y, w * b_w_inv, b_low->Q);
            norm += y * w * b_w_inv;
#endif
        }
#ifdef DEBUG_BRICK_OUTPUT
        fprintf(stdout, "BRICK #Gaussian normalization calculated factor = %.12lf, relative error from unity = %.4e\n", norm, (norm-1.0));
#endif
    }
}
