#include <math.h>
#include <assert.h>
#include <jibal_units.h>

#ifdef ERF_Q_FROM_GSL
#include <gsl/gsl_sf.h>
#endif
#include "brick.h"

extern inline double erfc_fast(double);

inline double erfc_fast(double x) {
    return exp(-1.0950081470333*x-0.75651138383854*x*x);
    /* Tsay, WJ., Huang, C.J., Fu, TT. et al. J Prod Anal 39, 259–269 (2013). https://doi.org/10.1007/s11123-012-0283-1 */
}

extern inline double erf_Q(double);

inline double erf_Q(double x) {
    return x < 0.0 ? 1.0-0.5*erfc_fast(-1.0*x/M_SQRT2) : 0.5*erfc_fast(x/M_SQRT2);
}

extern inline double erf_Q_new(double);

inline double erf_Q_new(double x) { /* Same as above, but without calling erfc_fast(). Saves a multiplication or two. */
    return x < 0.0 ? 1.0-0.5*exp(0.77428768622*x-0.37825569191*x*x) : 0.5*exp(-0.77428768622*x-0.37825569191*x*x);
}

#define ERFQ_A (-1.0950081470333/M_SQRT2)
#define ERFQ_B (-0.75651138383854*0.5)
#define ERFQ_CUTOFF 5.0 // 5.0 is enough

double erf_Q_optim(double x);

inline double erf_Q_optim(double x) { /* This turned out to be slower than expected */
    return x < (-ERFQ_CUTOFF) ? 1.0 : x > (ERFQ_CUTOFF) ? 0.0 : x < 0.0? 1.0-0.5*exp(-1.0*ERFQ_A*x+ERFQ_B*x*x) : 0.5*exp(ERFQ_A*x+ERFQ_B*x*x);
}

void erf_Q_test() {
    double x;
    int i;
    double sigma = 1.0;
    double x_L = 10.0;
    double x_H = 15.0;
    double sum = 0.0;
    double step = 0.03;
    for(i=0; i <10000; i++) {
        x = -50.0 + i*step;
#ifdef ERF_Q_FROM_GSL
        double fx = gsl_sf_erf_Q((x_H-x)/sigma);
        double fx_low = gsl_sf_erf_Q((x_L-x)/sigma);
#else
        double fx = erf_Q((x_H-x)/sigma);
        double fx_low = erf_Q((x_L-x)/sigma);
#endif
        double out = (fx_low-fx)/(x_H-x_L); /* One count, convolution of gaussian rectangle */
        sum += out*step; /* This is the integral! */
        fprintf(stdout, "%g %g %g %12.8lf\n", x, fx, out, sum);
    }
}

void bricks_calculate_sigma(const detector *det, const jibal_isotope *isotope, brick *bricks, size_t n_bricks) {
    for(size_t i = 0; i <= n_bricks; i++) {
        bricks[i].sigma = sqrt(bricks[i].S + detector_resolution(det, isotope, bricks[i].E) +  bricks[i].S_geo_x + bricks[i].S_geo_y);
    }
}

void brick_int2(gsl_histogram *h, const brick *bricks, size_t n_bricks, const double scale, const double sigmas_cutoff) {
    for(size_t i = 1; i <= n_bricks; i++) {
        const brick *b_high = &bricks[i-1];
        const brick *b_low = &bricks[i];
        //double E_diff_brick = b_high->E - b_low->E;
        //fprintf(stderr, "delta d = %.3lf, E_high = %.3lf, E_low = %.3lf), sigma_low = %.3lf, sigma_high = %.3lf, Q = %.3lf\n", (b_low->d - b_high->d)/C_TFU, b_high->E/C_KEV, b_low->E/C_KEV, sigma_low/C_KEV, sigma_high/C_KEV, b_low->Q);
        double E_cutoff_low = b_low->E - b_low->sigma * sigmas_cutoff; /* Low energy cutoff (brick) */
        double E_cutoff_high = b_high->E + b_high->sigma * sigmas_cutoff;
        for(size_t j = 0; j < h->n; j++) {
            if(h->range[j+1] < E_cutoff_low) /* High energy edge of histogram is below cutoff */
                continue;
            if(h->range[j] > E_cutoff_high) /* Low energy edge of histogram is above cutoff */
                break; /* Assumes histograms have increasing energy */
            const double E = (h->range[j] + h->range[j + 1]) / 2.0; /* Approximate gaussian at center bin */
            const double w = h->range[j + 1] - h->range[j];
            //const double y = (erf_Q((b_low->E - E) / b_low->sigma) - erf_Q((b_high->E - E) / b_high->sigma)) / (b_high->E - b_low->E);
            const double y = (erf_Q_new((b_low->E - E) / b_low->sigma) - erf_Q_new((b_high->E - E) / b_high->sigma)) / (b_high->E - b_low->E);
            h->bin[j] += scale * y * w * b_low->Q;
        }
    }
}

void brick_int(double sigma_low, double sigma_high, double E_low, double E_high, gsl_histogram *h, double Q) { /* Energy spectrum h, integrate over rectangular brick convoluted with a gaussian */
#ifdef OPTIMIZE_BRICK
    size_t lo;
    size_t hi;
    gsl_histogram_find(h, E_low, &lo);
    gsl_histogram_find(h, E_high, &hi);

#ifdef SMART
    double foo = (E_high-E_low)/(3.0); /* TODO: smart choice */
    int n = ceil(foo*5);

    if(n < 50)
        n = 50;
#else
    int n = 50; /* TODO: stupid choice */
#endif

    if(lo > n)
        lo -= n;
    else
        lo = 0;
    if(hi > h->n-n-1)
        hi = h->n;
    else
        hi += n;
    for(size_t i = lo; i <= hi; i++) {
#else
        for(size_t i = 0; i < h->n; i++) {
#endif
        double w = h->range[i+1] - h->range[i];
        double E = (h->range[i] + h->range[i+1])/2.0; /* Approximate gaussian at center bin */
#ifdef ERF_Q_FROM_GSL
        double y = (gsl_sf_erf_Q((E_low-E)/sigma)-gsl_sf_erf_Q((E_high-E)/sigma))/(E_high-E_low);
#else
#ifdef SIGMA_MEAN
        double sigma = (sigma_low+sigma_high)/2.0;
        double y = (erf_Q((E_low-E)/sigma)-erf_Q((E_high-E)/sigma))/(E_high-E_low);
#else
        double y = (erf_Q_optim((E_low-E)/sigma_low)-erf_Q_optim((E_high-E)/sigma_high))/(E_high-E_low);
#endif
#endif
        h->bin[i] += y*w*Q;
    }
}
