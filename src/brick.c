#include <math.h>
#include <assert.h>
#include <jibal_units.h>
#include <jibal_phys.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_math.h>
#include "brick.h"
#include "defaults.h"

#define BRICK_SQRT_2PI (2.50662827463100050242)

extern inline double erf_Q_fast(double x) { /* Approximative gaussian CDF */
    return x < 0.0 ? 1.0-0.5*exp(0.77428768622*x-0.37825569191*x*x) : 0.5*exp(-0.77428768622*x-0.37825569191*x*x);
    /* Based on: Tsay, WJ., Huang, C.J., Fu, TT. et al. J Prod Anal 39, 259–269 (2013). https://doi.org/10.1007/s11123-012-0283-1 */
}

void bricks_calculate_sigma(const detector *det, const jibal_isotope *isotope, brick *bricks, size_t last_brick) {
    for(size_t i = 0; i <= last_brick; i++) {
        bricks[i].S_det = detector_resolution(det, isotope, bricks[i].E);
        double S_gauss = bricks[i].S + bricks[i].S_det + bricks[i].S_geo_x + bricks[i].S_geo_y;
        bricks[i].S_sum = sqrt(S_gauss);
        for(int k = 0; k < MS_MIX_N; k++) {
            bricks[i].S_sum_ms[k] = bricks[i].ms_hwhm > 0.0 ? sqrt(S_gauss + bricks[i].ms.v[k]) : bricks[i].S_sum;
        }
    }
}

typedef struct {
    double w; /* Weight */
    double S_low, S_high; /* Widths (sigma) at low and high edges */
    int n_sub; /* Number of sub-bricks, zero if brick is not split */
    int degenerate; /* Brick has no width in energy, sub-bricks are gaussians */
    double S_inv_sub[BRICK_SPLIT_MAX]; /* Inverse of width of each sub-brick */
    double E_low, E_high; /* Energy range where this component (when split) contributes */
} brick_component;

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
        if(!b_low->valid) {
            continue;
        }
        const brick *b_high = &bricks[i-1];
        const int ms = (b_low->ms_hwhm > 0.0 || b_high->ms_hwhm > 0.0); /* Multiple scattering, convolution with a mixture of Gaussians */
        double ms_w[MS_MIX_N];
        double S_low = b_low->S_sum * sigmas_cutoff, S_high = b_high->S_sum * sigmas_cutoff; /* Cutoff distances */
        if(ms) {
            for(int k = 0; k < MS_MIX_N; k++) {
                ms_w[k] = (b_low->ms.w[k] + b_high->ms.w[k]) / 2.0; /* Components with same index are paired, both have variances in increasing order */
                if(ms_w[k] > 0.0) { /* Cutoff of a component depends on its weight, wide components in the tails have small weights */
                    double sigmas = GSL_MIN_DBL(sigmas_cutoff, sqrt(2.0 * log(GSL_MAX_DBL(ms_w[k] / BRICK_MS_COMPONENT_CUTOFF, 1.0))));
                    S_low = GSL_MAX_DBL(S_low, b_low->S_sum_ms[k] * sigmas);
                    S_high = GSL_MAX_DBL(S_high, b_high->S_sum_ms[k] * sigmas);
                }
            }
        }
        const int n_comp = ms ? MS_MIX_N : 1; /* Gaussian components of the energy spread */
        const double dE_brick = b_high->E - b_low->E;
        const double E_brick_min = GSL_MIN_DBL(b_low->E, b_high->E), E_brick_max = GSL_MAX_DBL(b_low->E, b_high->E);
        brick_component comp[MS_MIX_N];
        for(int k = 0; k < n_comp; k++) {
            brick_component *bc = &comp[k];
            bc->w = ms ? ms_w[k] : 1.0;
            if(bc->w <= 0.0) {
                continue;
            }
            const double S_l = ms ? b_low->S_sum_ms[k] : b_low->S_sum;
            const double S_h = ms ? b_high->S_sum_ms[k] : b_high->S_sum;
            const double S_min = GSL_MIN_DBL(S_l, S_h);
            bc->S_low = S_l;
            bc->S_high = S_h;
            bc->n_sub = 0; /* Zero: the brick is convoluted using different widths at the low and high edges */
            bc->degenerate = (fabs(dE_brick) < BRICK_DEGENERATE_WIDTH * S_min); /* Brick has (practically) no width, it is a gaussian */
            if(S_min > 0.0 && (bc->degenerate || sigmas_cutoff * fabs(S_h - S_l) > fabs(dE_brick))) {
                /* Brick is narrow compared to the change in width. Using different widths at the edges would result in large positive and negative lobes,
                 * so the brick is split into sub-bricks, each convoluted with a single width interpolated (in variance) to the middle of the sub-brick. */
                int n = (int)ceil(fabs(S_h - S_l) / (BRICK_SPLIT_SIGMA_STEP * S_min));
                bc->n_sub = GSL_MIN_INT(GSL_MAX_INT(n, 1), BRICK_SPLIT_MAX);
                for(int i_sub = 0; i_sub < bc->n_sub; i_sub++) {
                    const double t = (i_sub + 0.5) / bc->n_sub;
                    bc->S_inv_sub[i_sub] = 1.0 / sqrt(pow2(S_l) + (pow2(S_h) - pow2(S_l)) * t);
                }
                const double sigmas = ms ? GSL_MIN_DBL(sigmas_cutoff, sqrt(2.0 * log(GSL_MAX_DBL(bc->w / BRICK_MS_COMPONENT_CUTOFF, 1.0)))) : sigmas_cutoff;
                bc->E_low = E_brick_min - sigmas * GSL_MAX_DBL(S_l, S_h); /* Range where this component contributes */
                bc->E_high = E_brick_max + sigmas * GSL_MAX_DBL(S_l, S_h);
            }
        }
        double E_cutoff_low, E_cutoff_high; /* Low energy cutoff (gaussian), High energy cutoff (gaussian) */
        if(b_low->E < b_high->E) { /* Detected energy increasing as brick number increases */
            E_cutoff_low = b_low->E - S_low;
            E_cutoff_high = b_high->E + S_high;
        } else { /* Energy decreasing as brick number increases */
            E_cutoff_low = b_high->E - S_high;
            E_cutoff_high = b_low->E + S_low;
        }
        if(E_cutoff_high < emin) { /* Brick (before convolution) is entirely far below emin. Safe to ignore. */
            continue;
        }
        E_cutoff_low = GSL_MAX_DBL(E_cutoff_low, emin);
        const double b_w_inv = (dE_brick != 0.0) ? 1.0 / dE_brick : 0.0; /* inverse of brick width (in energy) */
        size_t ch_start = calibration_inverse(c, E_cutoff_low, h->n);
#ifdef DEBUG_BRICK_OUTPUT
        double norm = 0.0;
#endif
        for(size_t j = ch_start; j < h->n - 1; j++) {
            if(h->range[j] > E_cutoff_high) /* Low energy edge of histogram is above cutoff */
                break; /* Assumes histograms have increasing energy */
            const double E = (h->range[j] + h->range[j + 1]) / 2.0; /* Approximate gaussian at center bin */
            const double w = h->range[j + 1] - h->range[j];
            double dens = 0.0; /* Probability density (1/energy) at E */
            for(int k = 0; k < n_comp; k++) {
                const brick_component *bc = &comp[k];
                if(bc->w <= 0.0) {
                    continue;
                }
                double y;
                if(bc->n_sub == 0) {
                    y = (erf_Q((b_low->E - E) / bc->S_low) - erf_Q((b_high->E - E) / bc->S_high)) * b_w_inv;
                } else {
                    if(E < bc->E_low || E > bc->E_high) {
                        continue;
                    }
                    y = 0.0;
                    const double dE_sub = dE_brick / bc->n_sub;
                    for(int i_sub = 0; i_sub < bc->n_sub; i_sub++) {
                        const double S_inv = bc->S_inv_sub[i_sub];
                        const double E_sub_low = b_low->E + dE_sub * i_sub;
                        if(bc->degenerate) {
                            const double x = (E_sub_low + dE_sub / 2.0 - E) * S_inv;
                            y += exp(-0.5 * x * x) * S_inv / (BRICK_SQRT_2PI * bc->n_sub);
                        } else {
                            y += (erf_Q((E_sub_low - E) * S_inv) - erf_Q((E_sub_low + dE_sub - E) * S_inv)) * b_w_inv;
                        }
                    }
                }
                dens += bc->w * y;
            }
            double out = scale * dens * w * b_low->Q;
            h->bin[j] += out;
#ifdef DEBUG_BRICK_OUTPUT
            fprintf(stdout, "BRICK %3zu %4zu %9.3lf %9.3lf %12e %12e %12e %12e %12e\n", i, j, h->range[j] / C_KEV, h->range[j + 1] / C_KEV, out, h->bin[j], scale, dens * w, b_low->Q);
            norm += dens * w;
#endif
        }
#ifdef DEBUG_BRICK_OUTPUT
        fprintf(stdout, "BRICK #Gaussian normalization calculated factor = %.12lf, relative error from unity = %.4e\n", norm, (norm-1.0));
#endif
    }
}
